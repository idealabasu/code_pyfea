# -*- coding: utf-8 -*-
"""
Written by Daniel M. Aukes.
Email: danaukes<at>seas.harvard.edu.
Please see LICENSE.txt for full license.
"""
import sys
from cx_Freeze import setup, Executable
import os
from os.path import join,normpath,dirname
import shutil
import idealab_tools.setup_tools as st

# Remove the existing folders folder
shutil.rmtree("build", ignore_errors=True)
shutil.rmtree("dist", ignore_errors=True)

#import pydevtools.empty_project
import pyfea_examples
import importlib
#script = 'pydevtools.empty_project.empty_main_widget'
script = 'pyfea_design_tool.design_tool'
importlib.import_module(script)
script_path = sys.modules[script]
script2 = script_path.__file__

#program_base_dir = os.path.split(pydevtools.empty_project.__file__)[0]
#program_file = 'empty_main_widget.py'
#script2 = st.fix(program_base_dir,program_file)

package_data = {}
package_data['pyfea_design_tool'] = ['files/*']


packages = []
packages.append('pyfea')
packages.append('pyfea_examples')
packages.append('pyfea_examples.tube')
packages.append('pyfea_design_tool')

#packages.append('pydevtools.empty_project')
#packages.append("scipy")
#packages.append("scipy.sparse.csgraph._validation")
#packages.append('numpy')
#packages.append("scipy")
#packages.append('pyfea_examples')
#python_installed_directory = dirname(sys.executable)

zip_includes = []
include_files = []
includes = []
excludes = []

includes.append('pyfea')
includes.append('pyfea_examples')
includes.append('pyfea_examples.tube')
includes.append('pyfea_design_tool')
includes.append('OpenGL')
includes.append('meshio')
includes.append('pyqtgraph')
includes.append('pyqtgraph.debug')
includes.append('numpy.core._methods')
includes.append("scipy.sparse.csgraph._validation")
includes.append('pyqtgraph.ThreadsafeTimer')
includes.append('pyqtgraph.opengl.shaders')
#includes.append('mpl_toolkits')

#includes.append("scipy.spatial.ckdtree")
#includes.append("scipy.integrate.vode")
#includes.append("scipy.integrate.lsoda")

include_files.append((st.fix(st.python_installed_directory,'Library/bin/geos_c.dll'),'Library/bin/geos_c.dll'))
include_files.append((st.fix(st.python_installed_directory,'Library/bin/geos.dll'),'Library/bin/geos.dll'))
#necessary for qt5
include_files.extend(st.include_entire_directory(st.fix(st.python_installed_directory,'Library\\plugins\\platforms'),''))
#include_files.append(('C:\\bin\\gmsh\\gmsh.exe',''))
#include_files.append(('C:\\Users\\daukes\\scripts\\gmsh.exe',''))
include_files.append(('C:\\bin\\gmsh.exe',''))
include_files.extend(st.include_entire_directory(st.fix(st.python_installed_directory,'Library\\bin'),''))

#include_files.extend(st.include_entire_directory('C:\\Users\\daukes\\code\\pydevtools\\pydevtools\\empty_project\\files','files'))
#include_files.extend(include_entire_directory('C:\\Miniconda3\\tcl',''))

#zip_includes.extend(st.include_entire_directory('C:\\Users\\daukes\\code\\pydevtools\\pydevtools\\empty_project\\files','files'))
zip_includes.extend(st.include_entire_directory(st.fix(st.python_installed_directory,"Lib/site-packages/OpenGL"),"OpenGL"))
#zip_includes.append(('C:\\bin\\gmsh\\gmsh.exe',''))

excludes.append('gtk')
excludes.append('_gtkagg')
excludes.append('_tkagg')
excludes.append('bsddb')
excludes.append('curses')
#excludes.append('email')
excludes.append('pywin.debugger')
excludes.append('pywin.debugger.dbgcon')
excludes.append('pywin.dialogs')
excludes.append('tcl')
excludes.append('tk')
excludes.append('Tkconstants')
excludes.append('Tkinter')
excludes.append('babel')
#excludes.append('pandas')
excludes.append('notebook')
excludes.append('spyder')
excludes.append('ipython')
excludes.append('jupyter_client')
excludes.append('jupyter_core')


build_exe_options = {}
build_exe_options['packages']=packages
build_exe_options['includes']=includes
build_exe_options['excludes']=excludes
build_exe_options["include_files"]=include_files
build_exe_options["zip_includes"]=zip_includes
build_exe_options['include_msvcr']=True
#build_exe_options['icon']='logo_4_1_icon.ico'

bdist_msi_options = {}
uuid = '{6fcf34aa-ca63-45a5-8d9a-c59e7423f6d7}'
bdist_msi_options['upgrade_code']= uuid

setup_options = {}
setup_options['build_exe']=build_exe_options
setup_options['bdist_msi']=bdist_msi_options

base = "Win32GUI"


program_name = 'SoftDesignTool'

#"script":               #the name of the file containing the script which is to be frozen
#"initScript":           #the name of the initialization script that will be executed before the actual script is executed; this script is used to set up the environment for the executable; if a name is given without an absolute path the names of files in the initscripts subdirectory of the cx_Freeze package is searched
#"base":                 #the name of the base executable; if a name is given without an absolute path the names of files in the bases subdirectory of the cx_Freeze package is searched
#"path":                 #list of paths to search for modules
#"targetDir":            #the directory in which to place the target executable and any dependent files
#"targetName":           #the name of the target executable; the default value is the name of the script with the extension exchanged with the extension for the base executable
#"includes":             #list of names of modules to include
#"excludes":             #list of names of modules to exclude
#"packages":             #list of names of packages to include, including all of the package's submodules
#"replacePaths":         #Modify filenames attached to code objects, which appear in tracebacks. Pass a list of 2-tuples containing paths to search for and corresponding replacement values. A search for '*' will match the directory containing the entire package, leaving just the relative path to the module.
#"compress":             #boolean value indicating if the module bytecode should be compressed or not
#"copyDependentFiles":   #boolean value indicating if dependent files should be copied to the target directory or not
#"appendScriptToExe":    #boolean value indicating if the script module should be appended to the executable itself
#"appendScriptToLibrary":#boolean value indicating if the script module should be appended to the shared library zipfile
#"icon":                 #name of icon which should be included in the executable itself on Windows or placed in the target directory for other platforms
#"namespacePackages":    #list of packages to be treated as namespace packages (path is extended using pkgutil)
#"shortcutName":         #the name to give a shortcut for the executable when included in an MSI package
#"shortcutDir":          #the directory in which to place the shortcut when being installed by an MSI package; see the MSI Shortcut table documentation for more information on what values can be placed here.

exe = Executable(script = script2, 
                 base=base,
                 shortcutName=program_name,
                 shortcutDir="ProgramMenuFolder",
                 icon = 'python/pyfea_design_tool/files/logo_4_1_icon.ico', # if you want to use an icon file, specify the file name here
                 )


executables = []
executables.append(exe)

setup_arguments = {}
setup_arguments['name'] = program_name
setup_arguments['author'] = 'Dan Aukes'
setup_arguments['author_email'] = 'danaukes@asu.edu'
setup_arguments['version'] = '0.0.2'
setup_arguments['description'] = 'FEA-based design tool'
setup_arguments['executables'] = executables
setup_arguments['options'] = setup_options
setup_arguments['package_data'] = package_data

# patch to set environment variable for tcl and tk libraries
module = importlib.import_module('tcl')
p = list(module.__path__)[0]
os.environ['TCL_LIBRARY'] = p
os.environ['TK_LIBRARY'] = p

setup(**setup_arguments)        
