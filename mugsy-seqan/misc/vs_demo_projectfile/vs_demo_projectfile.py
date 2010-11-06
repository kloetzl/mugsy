import sys
import string
import os

BUILD_FLAGS_7 = '/I ".." /Op- /EHsc /D "DEBUG" /D "WIN32" /Zi /GR /W2 /Zc:wchar_t'
BUILD_FLAGS_8 = '/I ".." /EHsc /MTd /Zi /D "_CONSOLE" /D "_CRT_SECURE_NO_DEPRECATE" /D "DEBUG" /D "WIN32" /Zc:wchar_t /W2 /wd4996'
BUILD_FLAGS_9 = '/I ".." /EHsc /MTd /Zi /D "_CONSOLE" /D "_CRT_SECURE_NO_DEPRECATE" /D "DEBUG" /D "WIN32" /Zc:wchar_t /W2 /wd4996'

################################################################################

def addFile(name, build_flags):
    global vcproj1
    global vcproj2
    
    print ".",
    
    vcproj1 += '\t\t<Configuration Name="' + name + '|Win32" ConfigurationType="0">\n'
    vcproj1 += '\t\t\t<Tool Name="VCNMakeTool"\n\t\t\t\tBuildCommandLine="&quot;$(VCInstallDir)bin\\cl.exe&quot; &quot;$(ConfigurationName).cpp&quot; ' + build_flags + '" Output="&quot;$(ConfigurationName).exe&quot;"/>\n'
    vcproj1 += '\t\t</Configuration>\n'

    vcproj2 += '\t\t<File RelativePath=".\\' + name + '.cpp"></File>\n'
    

################################################################################

def scanDemos(search_path, build_flags):
    global vcproj1
    global vcproj2

    vcproj1 = ""
    vcproj2 = ""
    
    build_flags = build_flags.replace('"', '&quot;');
    
    for root, dirs, files in os.walk(search_path):
        if 'CVS' in dirs: dirs.remove('CVS')
        if '.svn' in dirs: dirs.remove('.svn')
        for file in files:
            pos = file.rfind(".")
            if pos < 0: continue
            if file[pos+1:] in ["cpp", "CPP"]:
                addFile(file[:pos], build_flags)

################################################################################

def createProject(search_path, version):
    global vcproj1
    global vcproj2

    f = open(version + ".vcproj")
    t = f.read()
    f.close()
    
    t = t.replace("####PLATZHALTER1####", vcproj1)
    t = t.replace("####PLATZHALTER2####", vcproj2)
    f = open(os.path.join(search_path, version + ".vcproj"), "w")
    f.write(t)
    f.close()
    

    f = open(version + ".sln")
    t = f.read()
    f.close()
    f = open(os.path.join(search_path, version + ".sln"), "w")
    f.write(t)
    f.close()

   
################################################################################
    

def main():
    global BUILD_FLAGS
    
    print "This script creates Visual Studio files (.vcproj and .sln) for SeqAn demos"
    print "The created files are stored in the demos folder"
    
    print "Visual Studio 7 (.net 2003) Files:",
    scanDemos("../../projects/library/demos", BUILD_FLAGS_7)
    createProject("../../projects/library/demos", "SeqAn_7")
    
    print
    print "Visual Studio 8 (.net 2005) Files:",
    scanDemos("../../projects/library/demos", BUILD_FLAGS_8)
    createProject("../../projects/library/demos", "SeqAn_8")
    
    print
    print "Visual Studio 9 (.net 2008) Files:",
    scanDemos("../../projects/library/demos", BUILD_FLAGS_9)
    createProject("../../projects/library/demos", "SeqAn_9")

    print
    print "Files created."
    
    
main()
