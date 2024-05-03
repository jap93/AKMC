#pragma once

#include <vector>

class WorkGroup
{

private:

    /** vector containing processors within the work group */
    std::vector<int> procs;

public:

    /**
     * @brief returns the cores within a WorkGroup
     * 
     * @return procs (std::vector<unsigned int> )
     */
    std::vector<int> getWorkGroup(void)
        {return procs;}

    /**
     * @brief Set the WorkGroup object to failed
     * 
     */
    void addCore(int core)
        {procs.push_back(core);}

    /**
     * @brief clear the WorkGroup
     * 
     */
    void clearWorkGroup(void)
        {procs.clear();}

    WorkGroup() {};
    virtual ~WorkGroup() {procs.clear();};

};
