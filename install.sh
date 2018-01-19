#!/bin/bash
# Set up environment ready for make3dgp
# usage:  ./install <install_dir>

if [ -z "$1" ]; then
	if [ -f ~/.bashrc ]; then
		echo "export PATH=\$PATH:${PWD}" >> ~/.bashrc
		echo "export GP3DSOURCELOC='${PWD}'" >> ~/.bashrc
	fi
	if [ -f ~/.bash_profile ]; then
		echo "export PATH=\$PATH:${PWD}" >> ~/.bashrc
		echo "export GP3DSOURCELOC='${PWD}'" >> ~/.bashrc
	fi
	if [ -f ~/.cshrc ]; then
		echo "set path = ( \$path ${PWD} )"  >> ~/.cshrc
		echo "setenv GP3DSOURCELOC '${PWD}'"  >> ~/.cshrc
	fi
else
	mkdir $1
	if [ $? == 0 ]; then
		fullpath=`readlink -f $1`
		cp *.f90 $fullpath
		cp make3dgp Makefile $fullpath
		if [ -f ~/.bashrc ]; then
			echo "export PATH=\$PATH:$fullpath" >> ~/.bashrc
			echo "export GP3DSOURCELOC='$fullpath'" >> ~/.bashrc
		fi
		if [ -f ~/.bash_profile ]; then
			echo "export PATH=\$PATH:$fullpath" >> ~/.bashrc
			echo "export GP3DSOURCELOC='$fullpath'" >> ~/.bashrc
		fi
		if [ -f ~/.cshrc ]; then
			echo "set path = ( \$path $fullpath )"  >> ~/.cshrc
			echo "setenv GP3DSOURCELOC '$fullpath'"  >> ~/.cshrc
		fi
	fi
fi
