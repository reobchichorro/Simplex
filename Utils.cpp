#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <fstream>
#include <math.h>
#include <filesystem>

// #define all(v) v.begin(),v.end()
// #define eps 1e-6

typedef std::string str;

// void printToFile(const str& filepath, const str& toPrint) {
//     std::ofstream file(filepath);
//     file << toPrint;
//     file.close();
// }

// void removeFileType(const str& filename, str& outname) {
//     outname = filename;
//     outname.erase(std::find(all(outname), '.'),outname.end());
// }

// void getFileType(const str& filename, str& filetype) {
//     filetype = filename.substr(filename.find('.'));
// }

// Timer::Timer(const std::string &msg): st(msg) {reset();}
// Timer::~Timer() {
//     if(st.size()!=0)
//         printNow(st);
// }
// void Timer::reset() { t0 = std::chrono::system_clock::now(); }
// void Timer::printNow(const std::string &st) const {  
//     std::cerr << "Time to " << st << " : "<< getTimeNow() << "\n";
// }
// double Timer::getTimeNow() const {
//     auto  t1 = std::chrono::system_clock::now();;
//     //clock_gettime(CLOCK_REALTIME, &t1);
//     return convertTimeMsecs(diff(t0,t1))/1000;
// }

// double Timer::convertTimeMsecs(const std::chrono::duration<double> td) const {        
//     return td.count()*1000;
// }

// std::chrono::duration<double> Timer::diff(TIME_T  start, TIME_T  end) const
// {
//     return end-start;
// }
