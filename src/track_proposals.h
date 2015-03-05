#ifndef TRACK_PROPOSALS_H
#define TRACK_PROPOSALS_H

#include <string>
#include <map>
#include <numeric>
#include <boost/circular_buffer.hpp>
#include <Rcpp.h>

#define ACCEPTTRACK 100

class track_acceptance {
    std::map<std::string, boost::circular_buffer<int>*> accept_circ_buff;    
    std::map<std::string, double> delta_map;

    public:
        void track_parameter(std::string param, double delta);
        double operator[](std::string param);
        void accept_reject(std::string param, int k);
        double accept_rate(std::string param);
        void erase_buffer();
        void modify_deltas();
}

#endif
