#! /bin/bash

ps -ef| grep fit_fg | grep -v grep | awk '{print "kill -9 " $2}' | sh

ps -ef| grep lmp_daily | grep -v grep | awk '{print "kill -9 " $2}' | sh
