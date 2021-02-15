// dkuzmy3 proj2

#include <iostream>
#include "Graph.h"

using namespace std;

int main(int argc, char* argv[]) {
    graph g;
    if(argc != 5) cout << "Usage: cpath <file> <source> <destination> <budget>" << endl;
    else {
        if(!g.read_file(argv[1]))
            std::cout << "could not open file '" << argv[1] << "'\n";
    }

    string source = argv[2];
    string destination = argv[3];
    int budget = atoi(argv[4]);

    cout << "Source: " << source << endl;
    cout << "Destination: " << destination << endl;
    cout << "Budget: " << budget << endl << endl;

    g.setSource(g.name2id(source));
    g.setDestination(g.name2id(destination));
    g.setBudget(budget);

    //g.display();

    g.dijsktraMod();
    g.printP();
    g.findPath();

    // report the cost/time curve p[i] = (cost i, time i) for all v

    // report the path from s to v, cost and time
    // can find all paths that fit the cost, choose the one with least time

    // design priority queue to hold <(total cost, total time), destination>


    return 0;
}
