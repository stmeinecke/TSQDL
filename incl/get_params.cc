// #include "get_params.hh"
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split

namespace  katana
{

double getCmdOption(char ** begin, char ** end, const std::string & option, double Default)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
	char* b;
	double d;
	d = strtod(*itr, &b);
	if (0 == d && *itr == b) 
		{	std::cout << "Input of Option ''" << option << "'' was wrong. Setting to default: " << option << " " << Default << std::endl;      // error handling.
   			return Default;
		}
	std::cout << "Set Option: "<< option << " = " << d << std::endl;
        return d;

    }
    return Default;
}

double getCmdOption(const std::string optionlist, const std::string &option, double Default, bool quiet)
{

std::vector<std::string> options;
std::vector<std::string>::iterator it;
double returnval=Default;
boost::split(options, optionlist, boost::is_any_of(", "), boost::token_compress_on);
//for(auto const &value: options)
//	{ std::cout << value << std::endl;}
it=std::find(options.begin(), options.end(), option); 
if(it!= options.end() and it+1!=options.end())
	{if(!quiet) std::cout<< "Set Option: " << option << " " << *(it+1) << std::endl;
	 returnval=stof(*(it+1));}
return returnval;
}

std::string getCmdOption(char ** begin, char ** end, const std::string & option, std::string Default)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
	std::cout << "Set Option: "<< option << " = " << *itr << std::endl;
        return *itr;

    }
    return Default;
}

std::string getCmdOption(const std::string optionlist, const std::string &option, std::string Default, bool quiet)
{
	std::vector<std::string> options;
	std::vector<std::string>::iterator it;
	std::string returnval=Default;
	boost::split(options, optionlist, boost::is_any_of(", "), boost::token_compress_on);
	it=std::find(options.begin(), options.end(), option); 
	if(it!= options.end() and it+1!=options.end())
		{if(!quiet) std::cout<< "Set Option: " << option << " " << *(it+1) << std::endl;
	 	returnval=*(it+1);}
return returnval;
}

bool getCmdOption_bool(char ** begin, char ** end, const std::string & option, bool Default)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end)
    {	bool Val=!Default;
	std::cout << "Set Option: "<< option << " = " << Val << std::endl;
        return Val;

    }
    return Default;
}

bool getCmdOption_bool(const std::string optionlist, const std::string &option, bool Default, bool quiet)
{
	std::vector<std::string> options;
	std::vector<std::string>::iterator it;
	bool returnval=Default;
	boost::split(options, optionlist, boost::is_any_of(", "), boost::token_compress_on);
	it=std::find(options.begin(), options.end(), option); 
	if(it!= options.end())
		{if(!quiet) std::cout<< "Set Option: " << option << " " << *(it+1) << std::endl;
	 	returnval=!Default;}
return returnval;
}





bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


//char* getCmdOption(char ** begin, char ** end, const std::string & option)
//{
//    char ** itr = std::find(begin, end, option);
//    if (itr != end && ++itr != end)
//    {
//        return *itr;
//    }
//    return 0;
//}
//
//bool cmdOptionExists(char** begin, char** end, const std::string& option)
//{
//    return std::find(begin, end, option) != end;
//}
}
