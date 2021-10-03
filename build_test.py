# build_test.py
# Copyright (c) 2013 Ahmed M. Hussein (amhussein4@gmail.com)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import os.path
import sys

common_sources = []
#preprocessor_sources = []
solver_sources = ["CASDevTest.cpp"]
#postprocessor_sources = []

compiler = "mpicxx"
include_arg = ""
link_arg = "libcasdev.a -lpng"

#preprocessor_output = ""
solver_output = "test"
#postprocessor_output = ""

object_extension = ".o"

print_only = False

def NeedToCompile(sourcefile,objectfile):
	# check if object file exists
	if(os.path.isfile(objectfile)):
		# check its time signature
		sourcetime = os.path.getmtime(sourcefile)
		objecttime = os.path.getmtime(objectfile)
		# leave a tolerance of 10 seconds to avoid any inaccurate time reporting
		if(sourcetime < (objecttime - 10)):
			return False
	return True

def Execute(command):
	if(print_only):
		print(command)
	else:
		print(command)
		os.system(command)

def GetObjectName(sourcename):
	objectname = sourcename[0:sourcename.rfind('.')] + object_extension
	return objectname

def CompileSource(sourcefile,debug):
	if(NeedToCompile(sourcefile,GetObjectName(sourcefile))):
		if(debug):	compile_string = compiler + "  -Wall -Wextra -pedantic -ansi -g " + include_arg + " -c " + sourcefile
		else:		compile_string = compiler + "  -Wall -Wextra -pedantic -ansi -O2 " + include_arg + " -c " + sourcefile
		Execute(compile_string)
		return True
	return False

def CompileSources(source_list,debug):
	for source in source_list:
		CompileSource(source,debug)

def LinkObjects(source_list,output):
	link_command = compiler + " "
	for source in source_list:
		link_command = link_command + GetObjectName(source) + " "
	link_command = link_command + " " + link_arg + " -o " + output
	Execute(link_command)

def Compile(debug):
	CompileSources(common_sources,debug)
	#CompileSources(preprocessor_sources,debug)
	CompileSources(solver_sources,debug)
	#CompileSources(postprocessor_sources,debug)

def Link():
	#LinkObjects(common_sources + preprocessor_sources,preprocessor_output)
	LinkObjects(common_sources + solver_sources,solver_output)
	#LinkObjects(common_sources + postprocessor_sources,postprocessor_output)

clean = False
debug = False

for i in range(1,len(sys.argv)):
	if(sys.argv[i] == "clean"):
		clean = True
	elif(sys.argv[i] == "debug"):
		debug = True

if(clean):
	Execute("rm *.o")
	Execute("rm " + solver_output)

Compile(debug)
Link()


