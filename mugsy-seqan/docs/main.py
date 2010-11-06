import dddoc
import dddoc_html
import sys
import os

# Command line: main.py <inpath> [ <module> ]*

if len(sys.argv) < 2: inpath = "../../projects/library"
else: inpath =  sys.argv[1]

print("Welcome to Dot.Dot.Doc")

buildfull = (len(sys.argv) < 3)
indexonly = not buildfull and (sys.argv[2] == 'indexonly')

if buildfull or indexonly:
    print("Scanning " + inpath + "...")
    os.path.normpath(inpath)
    dddoc.loadFiles(inpath)
else:
    i = 2;
    while (i < len(sys.argv)):
        modulepath = inpath + "/seqan/" + sys.argv[i]
        os.path.normpath(modulepath)
        print("Scanning " + modulepath + "...")
        dddoc.loadFiles(modulepath)
        i += 1


print("Scanning pages...")
dddoc.loadFiles("pages")

print("Scanning concepts...")
dddoc.loadFiles("concepts")


dddoc.DATA.init()

print("Create HTML Documentation...")
dddoc_html.createDocs("html", buildfull, indexonly)

if buildfull:
    print("Documentation created.")
else:
    print("Documentation updated.")

#raw_input("press return")