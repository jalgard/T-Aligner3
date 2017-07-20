#ifndef DNAUTIL_H
#define DNAUTIL_H

#include <string>
#include <locale>
#include <vector>
#include <sstream>
// for work with FASTA/FASTQ
#include <map>
#include <fstream>
#include <algorithm>

namespace libdna {

static char libdnaComplement(char a)
{
    switch(a) {
        case 'A' : return 'T';
        case 'T' : return 'A';
        case 'G' : return 'C';
        case 'C' : return 'G';
        case 'N' : return 'N';
        case 'a' : return 't';
        case 't' : return 'a';
        case 'g' : return 'c';
        case 'c' : return 'g';
        case 'n' : return 'n';
        default: return 'N';
    }
}

static std::string upDNA(const std::string& dna)
{
    std::string tmp=dna;
    for (size_t pos = 0; dna[pos] != '\0'; ++pos)
        tmp[pos] = std::toupper(tmp[pos]);
    return tmp;
}

static std::string lowDNA(const std::string& dna)
{
    std::string tmp=dna;
    for (size_t pos = 0; dna[pos] != '\0'; ++pos)
        tmp[pos] = std::tolower(tmp[pos]);
    return tmp;
}

static std::string comDNA(const std::string& dna)
{
    std::string tmp=dna;
    for(size_t pos = 0; dna[pos] != '\0'; ++pos)
	   tmp[pos]=libdnaComplement(dna[pos]);
    return tmp;
}

static std::string revDNA(const std::string& dna)
{
    std::string tmp="";
    for(int pos=dna.length()-1; pos >= 0; --pos)
        tmp += dna[pos];
    return tmp;
}

static std::string rcDNA(const std::string& dna)
{
    std::string tmp = revDNA(dna);
    for(size_t pos = 0; tmp[pos] != '\0'; ++pos)
	tmp[pos]=libdnaComplement(tmp[pos]);
    return tmp;
}


static std::vector<std::string>& libdnaSplit(const std::string& s, char delim, std::vector<std::string>& elems)
{
    std::stringstream ss(s); std::string item;
    while (std::getline(ss, item, delim)) elems.push_back(item);
    return elems;
}

static std::vector<std::string> Tokenize(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    libdnaSplit(s, delim, elems);
    return elems;
}

static std::string Token(const std::string& s, const std::string& start, const std::string& end)
{
    std::string token = "";
    std::size_t ss = s.find(start); std::size_t se = s.find(end, ss + start.size() + 1);
    if(ss != std::string::npos && se != std::string::npos && ss < se)
    {
	token = s.substr(ss + start.size(), se - ss - start.size());
    }
    return token;
}

static std::string MergeTokens(std::vector<std::string>& elems, size_t start, size_t end)
{
    std::string merged = "";
    if(end == 0) end = elems.size();
    for(size_t i = start; i < end; i++) merged+=elems[i];
    return merged;
}

static std::string i2s(const int& value)
{
    std::ostringstream tmp;
    tmp << value;
    return tmp.str();
}

static std::string d2s(const double& value, const int prec)
{
    std::ostringstream tmp;
    tmp.precision(prec); tmp  << value;
    return tmp.str();
}


// Generic FASTA reader in map
static void read_fasta(const char* file, std::map<std::string, std::string>& fae,
	bool name_space_split = false, bool upper = true)
{
    std::ifstream ifile(file); std::string line, last;
	while(std::getline(ifile, line, '\n'))
	{
		if(line[0] == '>')
		{
			last = line.substr(1); if(name_space_split) last = Tokenize(last,' ')[0];
		}
		else
        {
            if(upper) fae[last]+=upDNA(line);
            else fae[last]+=line;
        }
	}
}
// Generic FASTA reader in vector
static void read_fasta(const char* file, std::vector<std::vector<std::string> >& fae,
	bool name_space_split = false, bool upper = true)
{
    std::ifstream ifile(file); std::string line, last;
	while(std::getline(ifile, line, '\n'))
	{
		if(line[0] == '>')
		{
			last = line.substr(1); if(name_space_split) last = Tokenize(last,' ')[0];
            std::vector<std::string> entry;
            entry.push_back(last);
            entry.push_back("");
            fae.push_back(entry);
		}
		else
        {
            if(upper) fae.back()[1] += upDNA(line);
            else fae.back()[1] += line;
        }
	}
}

// Generic FASTQ reader, by default loads reads in map [read name -> read seq]
// can save memory, discarding read's name

static void read_fastq(const char* file, std::map<std::string, std::string>& fqe,
	bool rename_reads_to_ids = false)
{
    std::ifstream ifile(file); std::string line1, line2, line3, line4;
	unsigned int id_count=0;
	while(  std::getline(ifile, line1, '\n') &&
			std::getline(ifile, line2, '\n') &&
			std::getline(ifile, line3, '\n') &&
			std::getline(ifile, line4, '\n'))
	{
		if(line1[0] == '@' && line2.size() > 1 && line2.size() == line4.size())
		{
			if(!rename_reads_to_ids) fqe[line1] = line2;
			else fqe[i2s(++id_count)] = line2;
		}
	}
}

static std::string toDNA(const std::string& dna)
{
    std::string seq = "";
    for(size_t i = 0; i < dna.size(); i++)
    {
        char L = std::toupper(dna[i]);
        if(L == 'A' || L == 'G' || L == 'C' || L == 'T')
            seq += L;
    }
    return seq;
}

} // libdna

#endif
