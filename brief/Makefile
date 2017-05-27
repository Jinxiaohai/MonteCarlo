##########################################################
##                 Created  by  xiaohai                 ##
##                 Telphone : 18501781924               ##
##            E-mail : jinxiaohai@sinap.ac.cn           ##
##            E-mail : xiaohaijin@outlook.com           ##
##   Address : Shanghai Institute of Applied Physics    ##
##########################################################

EXE=brief.exe
OBJS=brief.o
CC=g++
CFLAGE= -Wall -g -std=c++11

$(EXE):$(OBJS)
	$(CC) $(OBJS) -o $(EXE)

$(OBJS):brief.cpp	;	$(CC) $(CFLAGE) -c brief.cpp

.PHONY:clean
clean:
	-rm $(OBJS) $(EXE)
