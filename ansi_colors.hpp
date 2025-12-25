#ifndef _ANSI_COLORS_HPP_
#define _ANSI_COLORS_HPP_

#include <string>

namespace ansi
{
	const std::string red     		= "\x1b[1;31m";
	const std::string green   		= "\x1b[1;32m";
	const std::string yellow  		= "\x1b[1;33m";
	const std::string blue    		= "\x1b[1;34m";
	const std::string magenta 		= "\x1b[1;35m";
	const std::string cyan    		= "\x1b[1;36m";
	const std::string white   		= "\x1b[1;37m";
	const std::string reset   		= "\x1b[0m";
	const std::string erase   		= "\x1b[2K\x1b\x0d";
	const std::string save_cursor   = "\x1b[s\x1b[1;36m";
	const std::string restore_cursor= "\x1b[u\x1b[0m";
}



#endif