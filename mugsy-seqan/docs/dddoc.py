import os
import copy
import string

################################################################################


class Line:
    def __init__(self, _nodes, _text):
        global ID
        
        self.nodes = _nodes
        self._text = _text.rstrip()
        self.id = ID
        ID += 1
        
    def __repr__(self):
        return "\nNode(" + str(self.nodes) + ",'" + self._text + "')"
        
    def name(self, index = 0):
        if len(self.nodes) > index:
            return self.nodes[index]
        else:
            return '(unknown)'
            
    def text(self):
        return self._text
        
            
class Data:
    def __init__ (self, _lines, _level):
        self.lines = _lines
        self.level = _level
        
    def __repr__(self):
        return "Data(" + str(self.lines) + ", " + str(self.level) + ")"

    def init(self):
        self.lines.sort(sortLineCompare)
        
        relations = self["globals.relations"].at_level(1).lines
        for relation in relations:
            to = relation.name(2)
            arr = splitName(relation.text())
            if len(arr) > 0:
                findRelation(self, arr, to)
                
        self.lines.sort(sortLineCompare)
        
    def __getitem__(self, str):
        return self.find(str)
        
    def find(self, str):
        arr = splitName(str)
        lines = []
        maxi = 0
        for line in self.lines:
            i = 0
            while (i < len(arr)) and (i + self.level < len(line.nodes)) and (arr[i] == line.nodes[i + self.level]):
                i += 1
            if i == len(arr):
                lines.append(line)
                
            elif maxi > i:
                break    
                
            maxi = i
            
        data = Data(lines, self.level + len(arr))
        return data

    def at_level(self, level = 0):
        lines = []
        for line in self.lines:
            if len(line.nodes) == self.level + level:
                lines.append(line)
            
        data = Data(lines, self.level + level)
        return data
        
    def sub_level(self, level = 1):
        lines = []
        for line in self.lines:
            if len(line.nodes) >= self.level + level:
                lines.append(line)
            
        data = Data(lines, self.level + level)
        return data
        
    def by_occ(self):
        lines = copy.copy(self.lines)
        lines.sort(sortLinesByOcc)
        data = Data(lines, self.level)
        return data        

    def empty(self):
        return len(self.lines) == 0
        
    def name(self, index = 0):
        if len(self.lines) > 0:
            return self.lines[0].name(index)
        return '(empty)'
        
    def text(self):
        str = ''
        for line in self.at_level(0).lines:
            if (str != ''): str += '\n'
            str += line.text()
        return str
        
    def keys(self, level = 0):
        dict = {}
        for line in self.lines:
            if len(line.nodes) > self.level + level:
                dict[line.nodes[self.level + level]] = 1
                
        arr = dict.keys()
        arr.sort()
        return arr
            
    def keys_by_occ(self, level = 0):
        dict = {}
        for line in self.lines:
            if len(line.nodes) > self.level + level:
                key = line.nodes[self.level + level]
                if not dict.has_key(key) or (dict[key] > line.id):
                    dict[key] = line.id

        dict2 = {}
        for key in dict:
            dict2[dict[key]] = key
            
        arr2 = dict2.keys()
        arr2.sort()
        
        arr = []
        for i in arr2:
            arr.append(dict2[i])
        
        return arr
        
         
################################################################################
         
def sortLineCompare(left, right):
    l = left.nodes
    r = right.nodes
    
    i = 0
    while (i < len(l)) and (i < len(r)):
        ret = cmp(l[i], r[i])
        if ret != 0: 
            return ret
        i += 1
        
    if len(l) < len(r): return -1
    elif len(l) > len(r): return 1
    elif left.id < right.id: return -1
    else: return 1
       
################################################################################
         
def sortLinesByOcc(left, right):
    if left.id < right.id: return -1
    else: return 1
         
################################################################################

def findRelation(data, arr, to):
    global DATA
    
    if len(arr) > 0:
        if (arr[0] == '*'): 
            sub_data = data.sub_level(1)
        else: 
            lines = []
            for line in data.lines:
                if line.name(data.level) == arr[0]:
                    lines.append(line)
            sub_data = Data(lines, data.level + 1)

        findRelation(sub_data, arr[1:], to)
        
    else:
        for line in data.at_level(0).lines:
            text = line.name(0) + '.' + line.name(1)
            entry = splitName(line.text())
            entry = entry[:2] + [to]
            DATA.lines.append(Line(entry, text))

################################################################################

DATA = Data([], 0)
ID = 0

################################################################################

def clearData():
    global DATA
    DATA = Data([], 0)
    global ID
    ID = 0

################################################################################

def loadDDDOCFile(filename):
    f = open(filename)
    text = f.readlines()
    f.close()

    return text

################################################################################

def loadCPPFile(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    
    ret = []
    
    #test for SEQAN_NO_DDDOC
    for line in lines:
        if line.find("SEQAN_NO_DDDOC") >= 0:
            return ret;
    
    
    incomment = False
    innextcomment = False
    inextract = False
    
    for line in lines:
        line = line.rstrip()
        str_line = ""
        if len(line) == 0:
            if not innextcomment and not incomment: 
                str_line = "."
            else: 
                str_line = " "
      
        while len(line) > 0 :
            if innextcomment:
                if line[len(line)-1] == "\\" :
                    if inextract: str_line += line[: len(line)-1]
                else:
                    if inextract: str_line += line
                    innextcomment = False
                break
                
            elif incomment:
                pos1 = line.find("*/")
                if pos1 < 0:
                    if inextract: str_line += line;
                    break;
                else:
                    if inextract: 
                        str_line += line[:pos1];
                        line = line[pos1 + 3:];
                    else:
                        line = line[pos1 + 2:];
                    incomment = False;
                    
            else:
                pos1 = line.find("/*")
                pos2 = line.find("//")
                pos3 = line.find('"')
                if (pos1 >= 0) and ((pos2 < 0) or (pos1 < pos2)) and ((pos3 < 0) or (pos1 < pos3)):
                    pos9 = line.find("*/", pos1 + 2)
                    if (len(line) > pos1 + 2):
                        inextract = (line[pos1 + 2] == "/") or (line[pos1 + 2] == "*")
                    else:
                        inextract = False
                    if pos9 < 0 : 
                        if inextract: str_line += line[pos1 + 3:]
                        incomment = True
                        break
                    else: 
                        if inextract: 
                            str_line += line[pos1 + 3: pos3]
                            line = line[pos9 + 3:]
                        else:
                            line = line[pos9 + 2:]
                        
                elif (pos2 >= 0) and ((pos3 < 0) or (pos2 < pos3)):
                    pos2b = pos2 + 2;
                    while ((pos2b < len(line)) and ((line[pos2b] == "/") or (line[pos2b] == "*"))):
                        pos2b += 1
                    inextract = (pos2b > pos2 + 2)
                    if line[len(line)-1] == "\\" :
                        if inextract: str_line += line[pos2b: len(line)-1]
                        innextcomment = True
                    else:
                        if inextract: str_line += line[pos2b:]
                    break
                    
                elif pos3 >= 0:
                    pos9 = line.find('"', pos3 + 2)
                    if pos9 < 0:
                        line = line[pos9+1:]
                        break
                    else:
                        break
                        
                else:
                    break

        ret = ret + [str_line]
    
    return ret

################################################################################

def testFileType(filename):
    pos = filename.rfind(".")
    if (pos >= 0): ext = filename[pos+1:]
    else: ext = ""

    if (ext in ["c", "C", "cpp", "CPP", "c++", "C++", "h", "H", "hpp", "HPP", "h++", "H++"]): return 2
    elif (ext in ["dddoc", "DDDOC"]): return 1
    else: return 0


################################################################################
    
def loadFile(filename):
    file_type = testFileType(filename)
    if (file_type == 2): return loadCPPFile(filename)
    elif (file_type == 1): return loadDDDOCFile(filename)
    else: raise "unknown file type"


################################################################################

def parseFile(filename):
    text = loadFile(filename)
    text.append('.')
    
    context = [[]]
    str = False
    
    for line in text:
        if line != '': 
            if line[0] == '.':
                parseString(str, context)
                str = line
            elif str:
                if str[len(str)-1] != '\n': str += '\n'
                str += line
                
################################################################################
                
def parseString(str, context):
    global DATA
    
    if not str or (str == '.'):
        return [[]]
        
    level = 0
    while (level < len(str)) and (str[level] == '.'): 
        level += 1
        
    str = str[level:]
        
    if (level < len(context)):
        del context[level:]
        
    if len(context) > 0:
        entry = copy.copy(context[len(context) - 1])
    else:
        entry = []

    key = ''
    text = ''
    
    pos = 0
    is_escaped = False
    c_quoted = ''
    while (pos < len(str)):
        c = str[pos]
        if c == "\x0d":
            pos += 1
            continue
        if c_quoted != "":
            if c_quoted == c: c_quoted = ""
            else: key += c                       
        elif is_escaped: 
            key += c
            is_escaped = False
        else:
            if c == '\\': is_escaped = True
            elif c in ['"', "'"]: c_quoted = c
            elif (c == ':'):
                key = str[0:pos]
                text = str[pos+1:]
                break
            else: key += c
                
        pos += 1    
        
    entry += splitName(key)
    DATA.lines.append(Line(entry, text))
    context.append(entry)


################################################################################

def splitName(line):
    pos = 0
    key = ""
    c_quoted = ""
    is_escaped = False
    li = []
    while (pos < len(line)):
        c = line[pos]
        if c_quoted != "":
            if c_quoted == c: c_quoted = ""
            else: key += c                       
        elif is_escaped: 
            key += c
            is_escaped = False
        else:
            if c == '\\': is_escaped = True
            elif c in ['"', "'"]: c_quoted = c
            elif c == '.':
                if key != "":
                    li.append(key)
                    key = ""
            elif c == '|':
                if key != "":
                    li.append(key)
                    key = ""
                rest = line[pos+1:]
                if len(rest)>0: li.append(rest)
                break;
            elif c != '\n':
                key += c
                
        pos += 1    
    
    if key != "": li.append(key)
    
    return li


################################################################################

def splitUrl(line):
    pos = 0
    key = ""
    c_quoted = ""
    is_escaped = False
    li = []
    while (pos < len(line)):
        c = line[pos]
        if c_quoted != "":
            if c_quoted == c: c_quoted = ""
            else: key += c                       
        elif is_escaped: 
            key += c
            is_escaped = False
        else:
            if c == '\\': is_escaped = True
            elif c in ['"', "'"]: c_quoted = c
            elif c == '|':
                if key != "":
                    li.append(key)
                    key = ""
            elif c != '\n':
                key += c
                
        pos += 1    
    
    if key != "": li.append(key)
    
    return li


################################################################################

def loadFiles(search_path):
    for root, dirs, files in os.walk(search_path):
        for file in files:
            path = os.path.join(root, file)
            if testFileType(path):
                parseFile(path)
        if 'CVS' in dirs:
            dirs.remove('CVS')

################################################################################
