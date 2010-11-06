import sys
import string
import os



################################################################################

def createProject(path, name, source_files, header_files, additional_files, version):

    f = open("template_" + version + ".vcproj")
    t = f.read()
    f.close()
    
    t = t.replace("####APPNAME####", name)
    t = t.replace("####HEADER_FILES####", header_files)
    t = t.replace("####SOURCE_FILES####", source_files)
    t = t.replace("####ADDITIONAL_FILES####", additional_files)
    f = open(os.path.join(path, name + "_" + version + ".vcproj"), "w")
    f.write(t)
    f.close()
    

    f = open("template_" + version + ".sln")
    t = f.read()

    t = t.replace("####APPNAME####", name)
    f = open(os.path.join(path, name + "_" + version + ".sln"), "w")
    f.write(t)
    f.close()

################################################################################

def createApp(name, path):
    
    source_files = ""
    header_files = ""
    additional_files = ""
    
    files = os.listdir(path)
    for file in files:
        if os.path.exists(os.path.join(path, file)) and not os.path.isdir(os.path.join(path, file)):
            pos = file.rfind(".")
            if pos < 0: ext = ""
            else: ext = file[pos+1:]
            if ext in ["cpp", "CPP"]:
                if file != "paramChooser.cpp":
                    source_files += "<File RelativePath=\".\\" + file + "\"></File>"
            else:
                if ext in ["h", "H", "hpp", "HPP"]:
                    header_files += "<File RelativePath=\".\\" + file + "\"></File>"
                else:
                    if not ext in ["vcproj", "sln", "ncb", "suo", "user"]:
                        additional_files += "<File RelativePath=\".\\" + file + "\"></File>"

    createProject(path, name, source_files, header_files, additional_files, "7")
    createProject(path, name, source_files, header_files, additional_files, "8")
    createProject(path, name, source_files, header_files, additional_files, "9")
    
    
################################################################################

def scanApps(search_path):
    dirs = os.listdir(search_path)
    if 'CVS' in dirs: dirs.remove('CVS')
    if '.svn' in dirs: dirs.remove('.svn')
    for file in dirs:
        dir = os.path.join(search_path, file)
        if os.path.isdir(dir):
            print dir
            createApp(file, dir)
    
   
################################################################################
    

def main():
    global BUILD_FLAGS
    
    print "This script creates Visual Studio files (.vcproj and .sln) for SeqAn apps"
    print "The created files are stored in each apps subfolder"
    print

    scanApps("../../projects/library/apps")
    
    print
    print "Files created."
    
    
main()