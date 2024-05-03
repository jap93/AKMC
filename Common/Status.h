#pragma once

#include "Constants.h"

class Status
{

private:

    unsigned int
        /** an integer indicating the status 0 = SUCCESS, 1 = FAILED, 3 = INTERRUPTED */
        simStatus = FAILED;

public:

    /**
     * @brief returns the status of the simulation
     * 
     * @return int 
     */
    int getStatus(void)
        {return simStatus;}

    /**
     * @brief Set the Status object to failed
     * 
     */
    void setStatusFailed(void)
        {simStatus = FAILED;}

    /**
     * @brief Set the Status object to success
     * 
     */
    void setStatusSuccess(void)
        {simStatus = SUCCESS;}

    /**
     * @brief Set the Status to indicate that it has been interruped, but will continue, by MPI
     * 
     */
    void setStatusInterrupt(void)
        {simStatus = INTERRUPTED;}

    /**
     * @brief Set the Status to indicate that it has been terminated by MPI
     * 
     */
    void setStatusTerminated(void)
        {simStatus = TERMINATED;}

    Status() {};
    virtual ~Status() {};

};