/**
 * @file logger.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-05-27
 *
 *
 * @section DESCRIPTION
 *
 * A poor man's logger class takes care of formatting and possibly
 * redirecting output to standard streams/files. It is a static class,
 * pretty much used as a name space, while keeping some stream handles
 * etc. in a unique private instance.
 *
 *
 *
 * @section LICENSE
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License and a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include "preProDebugFlags.h"

#include <iostream>
#include <fstream>
#include <string>
//#include <sstream>
#include <memory>
#include <list>
#include <regex>
#include <stdexcept>
#include <utility>//pair

using namespace std;

/**
 * To use this logger just include this header and then call one of
 * the static output member functions as described there.
 *
 */

class logger{

private:

  static unique_ptr<logger> log;

  typedef list<pair<ofstream*,string>> list_type;

  list_type myFiles;

  ostream* default_stream;


public:

  ostream* getOpenFile(string name=string("")){
    for(list_type::iterator it = myFiles.begin(); it != myFiles.end(); it++){
      if (it->first->good())
        if(name.compare("") || name.compare(it->second))
          return it->first;
    }
    return 0;
  }

  ostream* openNewFile(string name=string("")){
    if(name.compare(""))
      name="log.log";

    myFiles.push_back(make_pair(new ofstream(name),name));
    return myFiles.back().first;
  }


  static const regex removePath;

  logger():default_stream(&cout){}

  ~logger(){
    for(list_type::iterator it = myFiles.begin(); it != myFiles.end(); it++){
      delete it->first;
    }
    myFiles.clear();
  }


  /**
   *
   * Some use cases:
   *
   * just cout:
   * logger::write("Hello world!");
   *
   * My default usecases for debug logging:
   * logger::write("very important status msg", logLevel, __FILE__, logger::removePath);
   * logger::write("very important status msg", logLevel, __FILE__, logger::removePath, logger::toFile("my.log"));
   *
   * Use the logger to apply some weird formatting
   * logger::write("Why would you swap all 'o's and 'u's in your output?", regex("([ou])(.*)([ou])"), "$3$2$1" );
   *
   */
  static void write(string output, int logLevel=11, string header = string(), regex replace_in_header =regex(), string replace_by=string(), string deliminator=string(" : "), bool terminate_by_endl=true, ostream* toStream=0);
  void logWrite(string output, int logLevel=11, string header = string(), regex replace_in_header =regex(), string replace_by=string(), string deliminator=string(" : "), bool terminate_by_endl=true, ostream* toStream=0){

    if(logLevel > DEBUG){
      return;
    }
    if(toStream==0)
      toStream = default_stream;

    if(!toStream->good()){
      cerr << "logger switched back to std::cout due to stream error."<<endl;
    }
    toStream = &cout;

    //    cout << "formatting header '" << header << "' with regexp " << " by '" << replace_by << "'" << endl;
    string formattedHeader = regex_replace(header,replace_in_header,replace_by);
    //    cout << "result: " << formattedHeader << endl;

    (*toStream) << formattedHeader << deliminator << output;
    if(terminate_by_endl)
      (*toStream) << endl;

  }

  //some variants:

  /// no header version
  static void write(string output,int logLevel=11, regex replace = regex(), string replace_by = string(), bool terminate_by_endl=true, ostream *toStream=0);
  void logWrite(string output, int logLevel=11, regex replace = regex(), string replace_by = string(), bool terminate_by_endl=true, ostream *toStream=0){
    write("",logLevel,output,replace,replace_by,"",terminate_by_endl,toStream);
  }


  /// handle log files internally
  static ostream * toFile(string name = string());
  ostream * logToFile(string name = string()){
    ostream * re = log->getOpenFile(name);
    if (re==0)
      return log->openNewFile(name);
    return re;
  }


};


  const regex logger::removePath = regex("^.*/");

unique_ptr<logger> logger::log = unique_ptr<logger>(new logger());

void logger::write(string output, int logLevel, string header, regex replace_in_header, string replace_by, string deliminator, bool terminate_by_endl, ostream* toStream){
  logger::log->logWrite(output,logLevel,header,replace_in_header,replace_by,deliminator,terminate_by_endl,toStream);
}

void logger::write(string output, int logLevel, regex replace, string replace_by, bool terminate_by_endl, ostream* toStream){
  logger::log->logWrite(output,logLevel,replace,replace_by,terminate_by_endl,toStream);
}

ostream* logger::toFile(string name){
  return logger::log->logToFile(name);
}


#endif
