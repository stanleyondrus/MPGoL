#!/bin/bash

IS_VPN=$(ping -c 1 mastiff.cs.rpi.edu | grep icmp* | wc -l)

# SSH Credentials
USERNAME='rices'
SSH_ADDR='mastiff.cs.rpi.edu'

# Sync directories
REMOTE_DIR="ssh://$USERNAME@$SSH_ADDR/hw4"
LOCAL_DIR="/home/seanrice/repos/MPGoL/"

# Do we have VPN connection?
if [ $IS_VPN = '1' ]
then
    echo "Connected to mastiff over VPN"

    # Sync files
    echo "Syncing files"
    /usr/bin/unison -perms 0 -batch -auto -confirmbigdel=false $LOCAL_DIR $REMOTE_DIR
    echo

    # Compile
    echo "Compiling"
    ssh $USERNAME@$SSH_ADDR "make -C /home/parallel/spring-2019/rices/hw4"
    echo

    # Run
    echo "Run"
    ssh $USERNAME@$SSH_ADDR "/usr/local/mpich-3.2/bin/mpirun -np 4 hw4/assignment4-5"
    echo

    # Sync files
    echo "Syncing files"
    /usr/bin/unison -perms 0 -batch -auto -confirmbigdel=false $LOCAL_DIR $REMOTE_DIR
    echo
fi
