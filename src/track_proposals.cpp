#include "track_proposals.h"

void track_acceptance::track_parameter(std::string param, double delta)
{
    accept_circ_buff.insert(std::make_pair(param, 
                new boost::circular_buffer<int>(ACCEPTTRACK)));
    delta_map.insert(std::make_pair(param, delta));
}

double track_acceptance::operator[](std::string param)
{
    return delta_map[param];
}

void track_acceptance::accept_reject(std::string param, int k)
{
    accept_circ_buff[param]->push_back(k);
}

double track_acceptance::accept_rate(std::string param)
{
    double acceptances = std::accumulate(accept_circ_buff[param]->begin(),
                                         accept_circ_buff[param]->end(), 0);
    double rate = acceptances / ACCEPTTRACK;

    return rate;
}

void track_acceptance::erase_buffer()
{
    for (auto itr=accept_circ_buff.begin(); 
         itr!=accept_circ_buff.end(); ++itr) {
        delete itr->second;
        accept_circ_buff.erase(itr);
    }
}

void track_acceptance::modify_deltas()
{
    for (auto itr=delta_map.begin(); itr!=delta_map.end(); ++itr) {
        double rate = accept_rate(itr->first);
        double current_delta = delta_map[itr->first];
        double delta_change;

        if (rate > 0.5) {
            // too many acceptances, increase size of delta
            delta_change = 1 + rate - 0.5;
            delta_map[itr->first] = current_delta * delta_change;
        } else if (rate < 0.15) {
            // too few acceptances, decrease size of delta
            delta_change = 1 - (0.15 - rate);
            delta_map[itr->first] = current_delta * delta_change;
        }
    }
}
