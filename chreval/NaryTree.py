#!/usr/bin/env python3

"""
Main Module:   Threads.py 

Objects:       Node
               NodeAnalysis
               TreeBuilder

Author:        Wayne Dawson
creation date: 190309
last update:   200211 (upgraded to python3), 191016
version:       0.0



an example of an N-ary tree in python

I was helped by the following sources on initially gettings started with this 

sources:

1. stackoverflow:
https://stackoverflow.com/questions/46570249/python-n-ary-tree-for-minimax
https://stackoverflow.com/questions/13730404/traversing-a-non-binary-tree


2. Problem Solving in Data Structures & Algorithms Using Python: Programming Interview Guide
by Hemant Jain

3. https://ruslanspivak.com/lsbasi-part1/ to https://ruslanspivak.com/lsbasi-part14/

"""

import sys
#from types import ListType
"""@

With python3, must use int, float, list, tuple; not the functions
intType, floatType, ListType, TupleType, etc, """
# Tests that can be run to understand this package.

TEST1 = True # False # 
TEST2 = False # True # 




class Node(object):
    def __init__(self, name=None, children=None):
        self.name = name
        self.children = children or [] # if None, assignes []
    #
    
    def __str__(self):
        s = "name(%s), " % self.name
        for ch in self.children:
            if not ch == None:
                s += "%s " % ch
            else:
                s += "() "
            #
        #
        return s
    #
    
    def __repr__(self):
        return self.__str__()
    #
    
                     
    
    
    # this is the proposed solution for getting the height of the
    # tree.
    
    def height(self):
        if not self.children:  # base case
            return 1
        else:                  # recursive case
            return 1 + max(child.height() for child in self.children)
        #
    #
    
    def disp_levels(self, level):
        s = ''
        s += "level %d: %s" % (level, self.name)
        for child in self.children:
            s += child.disp_levels(level+1)
        #
        return s
    #
    
    # Not particularly useful, frankly, but anyway, it adds the same
    # thing to all the children in Node class. Hence, it is not
    # particularly useful in general.
    def addToAllChildren(self, name):
        if self.children:
            # print ("addToAllChildren: children: ", self.children)
            for node in self.children: # loop over all children of root
                # print ("addToAllChildren(in loop): ", name)
                node.addToAllChildren(name)         # add all the new names
            #
        else:
            # print ("extend", name)
            self.name.extend(name)
        #
    #
    
    # this is a generator, but I don't really like it very much
    def preorder(self):
        yield self.name # yield is used as a generator
        for child in self.children:
            for descendent in child.preorder():
                yield descendent
            #
            
            # yield from child.preorder() # Python 3.3 only!
        #|endfor
        
    #
    
#

class NodeAnalysis(object):
    def __init__(self):
        self.dvlevels = {}
        self.dstruct  = {}
        self.max_vlev = 0
        self.max_hpos = 0
        self.max_height = 0
        self.debug_structLayout = False # True # 
        self.debug_NodeAnalysis = False # True # 
    #
    
    def vertical_levels(self, vlevel, tree):
        self.max_height = tree.height()
        self.max_vlev = 0
        self.max_hpos = 0
        self.dvlevels = {}
        self.add_vlevel(vlevel, tree.name)
        # print ("tree.name:     ", tree.name)
        # print ("tree.children: ", tree.children)
        for child in tree.children:
            self.vertical_levels(vlevel+1, child)
        #
    #
    
    def add_vlevel(self, vlevel, data):
        if vlevel in self.dvlevels:
            # print ("1 vlevel %d, %s" % (vlevel, data))
            self.dvlevels[vlevel] += [data]
        else:
            # print ("0 vlevel %d, %s" % (vlevel, data))
            self.dvlevels.update({vlevel : [data]})
            if self.max_vlev < vlevel:
                self.max_vlev = vlevel
            #
        #
    #
    
    def disp_vlevels(self):
        print (self.dvlevels)
        for i in range(0, self.max_vlev + 1):
            s = "%2d:   " % i
            for vleveli in self.dvlevels[i]:
                for hdata in vleveli:
                    s += "%s, " % hdata
                #
            #
            print (s)
        #
    #
    
    def structLayout(self, tree):
        set_arr = [0]
        self.max_height = tree.height()
        self.dstruct = {}
        if self.debug_structLayout:
            print ("Enter structLayout: ")
            print ("maximum height = %d, level = %d" % (self.max_height, set_arr[0]))
        #
        
        if tree.name == None:
            if self.debug_structLayout:
                print ("assign empty root: ")
            #
            self.dstruct.update({ 0 : [[]] })
        else:
            self.dstruct.update({ 0 : [[tree.name]] })
            nlen = len(tree.children)
            if self.debug_structLayout:
                print ("nlen = %d" % nlen)
                print ("tree.children -> %s" % tree.children)
            #
            if nlen == 0:
                if self.debug_structLayout:
                    print ("assign empty children: ")
                #
                self.dstruct.update({ 1 : [[]] })
            else:
                level   = 1
                new_arr = [level, 0]
                
                v = []
                for pos in range(0, nlen):
                    new_arr[level] = pos
                    v += [tree.children[pos].name]
                    if self.debug_structLayout:
                        print ("new_arr(level = %d, pos = %d): " % (level, pos), new_arr)
                    #
                    self.structLayoutUtil(level, pos, new_arr, tree.children[pos])
                #
                if level in self.dstruct:
                    self.dstruct[level] += [v]
                else:
                    self.dstruct.update({ level : [v] })
                #
                if self.debug_structLayout:
                    print (self.dstruct)
                #
            #
        #
    #
    
    def structLayoutUtil(self, level, pos, set_arr, node):
        if self.debug_structLayout:
            print ("Enter structLayoutUtil: ")
            print ("set_arr (the goal): ", set_arr)
            print ("current level:      ", level)
            print ("current node name:  ", node.name)
            print ("number of children: ", len(node.children))
        #
        if set_arr[0] > self.max_height:
            print ("ERROR: maximum height of the tree (%d) is exceeded" \
                % (self.max_height - 1))
            sys.exit(0)
        #
        nlen = len(node.children)
        if nlen == 0: # reached a leaf
            level += 1
            if self.debug_structLayout:
                print ("n = 0: level = %d, len(set_arr) = %d, set_arr = %s" \
                    % (level, len(set_arr), set_arr))
                print ("assign empty children: ")
            #
            if level < self.max_height:
                # print ("1 assign empty root: ")
                for vlev in range(level, self.max_height):
                    if vlev in self.dstruct:
                        self.dstruct[vlev] += [[]]
                    else:
                        self.dstruct.update({ vlev : [[]] })
                    #
                #
            #
            if self.debug_structLayout:
                print (self.dstruct)
            #
        else:
            # print (self.dstruct)
            new_arr = []
            for n in set_arr:
                new_arr += [n]
            #
            
            if self.debug_structLayout:
                print ("n > 0: level = %d, len(new_arr) = %d, %s" \
                    % (level, len(set_arr), new_arr))
            #
            new_arr += [0]
            level += 1
            new_arr[0] = level
            
            v = []
            for pos in range(0, nlen):
                new_arr[level] = pos
                if self.debug_structLayout:
                    print ("       level = %d, len(new_arr) = %d, %s" \
                        % (level, len(new_arr), new_arr))
                #
                v += [node.children[pos].name]
                self.structLayoutUtil(level, pos, new_arr, node.children[pos])
            #
            if level in self.dstruct:
                self.dstruct[level] += [v]
            else:
                self.dstruct.update({ level : [v] })
            #
            # print (self.dstruct)
        #
    #
    
    def disp_structLayout(self):
        if self.debug_structLayout:
            print ("Enter disp_struct: ")
            print ("dstruct:        ", self.dstruct)
            print ("maximum height: ", self.max_height)
        #
        
        s = ''
        for vlev in range(0, self.max_height):
            s += "%2d:   " % vlev
            n = len(self.dstruct[vlev])-1
            for hpos in range(0, n):
                s += "%s | " % (self.dstruct[vlev][hpos])
            #|endfor
            
            if vlev == self.max_height - 1:
                s += "%s" % (self.dstruct[vlev][n])
            else:
                s += "%s\n" % (self.dstruct[vlev][n])
            #
            
        #|endfor
        
        return s
    #
#



class TreeBuilder(object):
    def __init__(self, name = None, children = None):
        self.tree = Node(name, children)
    #

    def howDeep(self, set_arr):
        depth, rslt = self.howDeepUtil(0, 0, set_arr)
        return depth
    #
    
    def howDeepUtil(self, itr, fwrd, set_arr):
        debug_howDeepUtil = False # True # 
        if itr > 40:
            # obviously something must be wrong if it ends up here!
            print ("ERROR(howDeepUtil): too many levels")
            sys.exit(0)
        #
        rslt = self.existsNode(set_arr)
        depth = set_arr[0]
        new_depth = depth
        if rslt:
            # have to permanently set the direction, ... otherwise, we
            # will osciliate between set_arr++ and set_arr--.
            use_plus = False
            if fwrd >= 0:
                fwrd = 1
                use_plus = True
            #
            if debug_howDeepUtil:
                print ("rslt, fwrd = ", rslt, fwrd)
            #
            if use_plus:
                # make a hard copy if set_arr and go deeper along the
                # chain (in essence, set_arr++)
                chk_arr = []
                for i in range(0, len(set_arr)):
                    chk_arr += [set_arr[i]]
                #
                chk_arr += [0]
                chk_arr[0] += 1
                if debug_howDeepUtil:
                    print ("+     set_arr = ", set_arr, itr)
                    print ("+ new chk_arr = ", chk_arr)
                #
                
                # search deeper
                new_depth, rsltd = self.howDeepUtil(itr + 1, fwrd, chk_arr)
                if rsltd:
                    depth = new_depth
                else:
                    new_depth = depth
                #
                if debug_howDeepUtil:
                    print ("new_depth(%d), depth(%d)" % (new_depth, depth))
                    if new_depth == depth:
                        print ("+ no change")
                    #
                #
            #
        
        else:
            # have to permanently set the direction, ... otherwise, we
            # will osciliate between set_arr++ and set_arr--.
            use_minus = False
            if fwrd <= 0:
                fwrd = -1
                use_minus = True
            #
            if debug_howDeepUtil:
                print ("rslt, fwrd = ", rslt, fwrd)
            #
            if use_minus:
                # make a hard copy if set_arr and return up along the
                # chain (in essence, set_arr--)
                chk_arr = []
                for i in range(0, len(set_arr)):
                    chk_arr += [set_arr[i]]
                #
                chk_arr[0] -= 1
                n = len(chk_arr)-1
                del chk_arr[n]
                if debug_howDeepUtil:
                    print ("      set_arr = ", set_arr, itr)
                    print ("- new chk_arr = ", chk_arr)
                #
                new_depth, rsltm = self.howDeepUtil(itr + 1, fwrd, chk_arr)
                if debug_howDeepUtil:
                    print ("new_depth(%d), depth(%d)" % (new_depth, depth))
                    if new_depth == depth:
                        print ("- no change")
                    #
                #
            #
        #
        return new_depth, rslt
    #

    def copyTree(self):
        tree2 = TreeBuilder()
        tree2.tree = self.copyTreeUtil(self.tree)
        return tree2
    #
    
    def copyTreeUtil(self, curr):
        if curr != None:
            temp = Node(curr.name)
            for k in range(0, len(curr.children)):
                temp.children += [self.copyTreeUtil(curr.children[k])]
            #
            return temp
        else:
            return None
        #
    #
    
    def numNodes(self):
        if self.tree.name == None:
            return 0
        else:
            return 1 + self.numNodesUtil(self.tree)
        #
    #
    
    def numNodesUtil(self, curr):
        if curr == None:
            return 0
        
        else:
            n = 0
            for k in range(0, len(curr.children)):
                n += 1 + self.numNodesUtil(curr.children[k])
            #
            
            return n
        #
    #
    
    def numLeafNodes(self):
        return self.numLeafNodesUtil(self.tree)
    #
    
    def numLeafNodesUtil(self, curr):
        if curr == None:
            return 0
        
        if len(curr.children) == 0:
            return 1
        
        else:
            n = 0
            for k in range(0, len(curr.children)):
                n +=  self.numLeafNodesUtil(curr.children[k])
            #
            
            return n
        #
    #
    
    def deleteChild(self, ch_arr):
        print ("Enter deleteChild: ", ch_arr)
        # print (self.tree)
        self.tree = self.deleteChildUtil(self.tree, 1, ch_arr)
    #
    
    def deleteChildUtil(self, node, level, ch_arr):
        # print ("Enter deleteChild: ", ch_arr)
        # print (node)
        if node != None:
            if level == ch_arr[0]:
                print ("level %d reached: %s, remove %d (last position)" % (ch_arr[0], ch_arr, ch_arr[level]))
                if node.children == None:
                    return None
                else:
                    del node.children[ch_arr[level]]
                    # print ("returning: ", node)
                    return node
                #
            #
            else:
                print ("level %d interum: %s" % (ch_arr[0], ch_arr[level]))
                node.children[ch_arr[level]]  = self.deleteChildUtil(node.children[ch_arr[level]], level+1, ch_arr)
            #
        #
        # print ("exiting:", node)
        return node
    #
    
    
    
    def insertChildren(self, ch_arr, children):
        self.tree = self.insertChildrenUtil(self.tree, 1, ch_arr, children)
    #
    
    def insertChildrenUtil(self, node, level, ch_arr, children):
        if level == ch_arr[0]:
            print ("level reached: ", ch_arr, ch_arr[level])
            # print (node)
            node.children = children
            # print (node)
            
        else:
            print ("level interum: ", ch_arr, ch_arr[level])
            node.children[ch_arr[level]]  \
                = self.insertChildrenUtil(node.children[ch_arr[level]],
                                          level+1,
                                          ch_arr,
                                          children)
        #
        
        return node
    #
        
    def add_newNode(self, set_arr, children):
        depth = self.howDeep(set_arr)
        n = len(set_arr) - 1
        if depth == set_arr[0]:
            print ("ERROR: slot %s is already filled!" % (set_arr))
            print ("     + chose an empty mode or directly")
            print ("       append data at this position")
            depth = -1
            
        elif depth < set_arr[0] - 1:
            print ("ERROR: slot %s is too deep (no neighoring node!" \
                % (set_arr))
            print ("       depth of node list along %s, %d, requested %d" \
                % (set_arr, depth, set_arr[0]))
            print ("     + fill in the empty modes between %d and %d first" \
                % (depth, set_arr[0] - 1))
            depth = -2
            
        elif set_arr[n] != 0:
            print ("ERROR: slot %s -> last index (%d) should be  0!" % (set_arr, set_arr[n]))
            print ("       suggest modifying string accordingly")
            depth = -3
            
        else:
            
            print ("requested insert at %s granted" % (set_arr))
            self.insertChildren(set_arr, children)
        #
        
        return depth
    #
        
    
    def check_coordinates(self, set_arr):
        flag_pass = True
        # print (set_arr[0], (len(set_arr)-1))
        if not set_arr[0] == (len(set_arr)-1):
            sa_len = len(set_arr)
            sa_set = set_arr[0]
            if (len(set_arr) > set_arr[0] + 1):
                print ("Error: too many entries to build tree array: ", set_arr)
                print ("       array should have %d elements, but has %d" \
                    % (sa_set + 1, sa_len))
                flag_pass = False
                
            elif (len(set_arr) < set_arr[0] + 1):
                print ("Error: insufficient information to build tree array: ", set_arr)
                print ("       array should have %d elements, but has only %d" \
                    % (sa_set + 1, sa_len))
                flag_pass = False
                
            else:
                print ("ERROR: unknow reason but arrangement of input coordinates is wrong")
                print ("                   set_arr = ", set_arr)
                print ("            len of set_arr = ", len(set_arr))
                print ("requested depth of set_arr = ", set_arr[0])
                flag_pass = False
            #
            
        #
        
        return flag_pass
    #
    
    
    # This sort of stuff was a pain to build, but it does inform if
    # the node exists and is assigned.
    def existsNode(self, set_arr):
        vlev = self.tree.height()
        # print ("maximum height: ", vlev, set_arr[0])
        if set_arr[0] >= vlev and vlev > 0:
            # print ("level too deep")
            return False
        #
        
        if not self.check_coordinates(set_arr):
            sys.exit(0)
        #
        
        if set_arr[0] == 0:
            # print ("0")
            if self.tree.name == None:
                return False
            else:
                return True
            #
            
        #
        
        nlen = len(self.tree.children)
        # print (nlen, self.tree.children)
        if nlen == 0:
            return False
        
        else:
            level = 1
            pos   = set_arr[level]
            # print ("new level = %d, nlen = %d, pos = %d" % (level, nlen, pos))
            if nlen > pos:
                return self.existsNodeUtil(set_arr, level, self.tree.children[pos])
            else:
                return False
            #
            
        #
        
    #
    
    def existsNodeUtil(self, set_arr, level, node):
        # print ("Enter existsNodeUtil: ")
        # print ("set_arr (the goal): ", set_arr)
        # print ("current level:      ", level)
        # print ("current node name:  ", node.name)
        # print ("number of children: ", len(node.children))
        
        # first, I just look to see if I can get some height
        if set_arr[0] == level:
            pos = set_arr[level]
            nlen = len(node.children)
            # print ("1 pos: %d, level = %d" % (pos, level))
            return True
        
        else:
            if set_arr[0] > level:
                nlen = len(node.children)
                # print (nlen, node.children)
                level += 1
                pos = set_arr[level]
                # print ("2 pos: %d, level = %d" % (pos, level))
                
                if nlen == 0:
                    return False
                elif nlen > pos:
                    return self.existsNodeUtil(set_arr, level, node.children[pos])
                else:
                    return False
                #
                
            else:
                print ("ERROR!!!!!!!!!!! existsNodeUtil")
                sys.exit(1)
            #
            
        #
    #
    
    
    def readNode(self, set_arr):
        vlev = self.tree.height() 
        # print ("maximum height: ", vlev, set_arr[0])
        if set_arr[0] >= vlev and vlev > 0:
            # print ("level too deep")
            return None
        #
        
        if not self.check_coordinates(set_arr):
            sys.exit(0)
        #
        
        if set_arr[0] == 0:
            # print ("0")
            if self.tree.name == None:
                return None
            else:
                return self.tree.name
            #
            
        #
        
        nlen = len(self.tree.children)
        # print (nlen, self.tree.children)
        if nlen == 0:
            return None
        
        else:
            level = 1
            pos   = set_arr[level]
            # print ("new level = %d, nlen = %d, pos = %d" % (level, nlen, pos))
            if nlen > pos:
                return self.readNodeUtil(set_arr, level, self.tree.children[pos])
            else:
                return None
            #
            
        #
        
    #
    
    def readNodeUtil(self, set_arr, level, node):
        # print ("Enter readNodeUtil: ")
        # print ("set_arr (the goal): ", set_arr)
        # print ("current level:      ", level)
        # print ("current node name:  ", node.name)
        # print ("number of children: ", len(node.children))
        
        # first, I just look to see if I can get some height
        if set_arr[0] == level:
            pos = set_arr[level]
            nlen = len(node.children)
            # print ("1 pos: %d, level = %d" % (pos, level))
            return node.name
        
        else:
            if set_arr[0] > level:
                nlen = len(node.children)
                # print (nlen, node.children)
                level += 1
                pos = set_arr[level]
                # print ("2 pos: %d, level = %d" % (pos, level))
                
                if nlen == 0:
                    return None
                elif nlen > pos:
                    return self.readNodeUtil(set_arr, level, node.children[pos])
                else:
                    return None
                #
                
            else:
                print ("ERROR!!!!!!!!!!! readNodeUtil")
                sys.exit(1)
            #
            
        #
    #


    
#    


"""
These are here because I was trying to figure out how to identify
if a variable was defined. This, I was thinking, would be useful for
determining if the variable was defined. However, as you see from my
comments, it looks like the only way to do this kind of test is do it
locally. Anyway, I found other ways to accomplish this, but it seems
useful to keep this information around.

# information for these tests (somewhat helpful) came from this source
# https://stackoverflow.com/questions/11556234/how-to-check-if-a-list-exists-in-python/11559451

"""
def testIsList(ima_list):
    # from types import ListType # defined above
    flag_IsList = False
    

    """
    There is a trick that can work if we could get the name of the
    variable to be expressed.

    mylist=[1,2,3]
    'mylist' in locals().keys()
    
    If 'mylist' is not in locals(), then it will deliver
    False. However, the difficulty is in defining this string. It
    seems like it can only be used locally for a fixed variable. The
    other thing is that a call to this function probably already
    produces an error because the program looks up locals() and sees
    that it is not there. So a call to this function already kills
    things if the variable doesn't exist.
    
    It may be the only thing that will work is to do the following
    
       try:
          #test the variable
          print ("defined")
       except NameError:
          print ("undefined")
       #
    
    Here I test this second idea, the first idea works as I explain
    above.

    """
    flag_IsList = False
    print ("first try:")
    try:
        # not myList is not defined yet!!!
        if type(myList) is list:
            flag_IsList = True
        #
    except NameError:
        print ("myList is undefined")
    #
    
    """
    apparently, this works, but it has to be inline with the variable
    being checked. At any rate, when we pass a variable to a function,
    the function looks up in its tables what is there, and if it is
    not found, you will get an error (for a list in a tree.children
    when the children's children don't exist). So it still require
    in-place.

    """    
    
    try:
        # for an existing/nonexisting list
        if type(ima_list) is list:  # was "if type(ima_list) is ListType:"
            flag_IsList = True
            print ("Attribute exists")
        #
        
    except (NameError, IndexError, AttributeError) as error:
        print (error)
        print ("problems")
        flag_IsList = False
    #
    return flag_IsList
#

def testIsTuple(ima_tuple):
    flag_isTuple = False
    # from types import TupleType # defined above
    try:
        # for an existing/nonexisting tuple
        if type(ima_tuple) is tuple: 
            flag_isTuple = True
            print ("Attribute exists")
        #
        
    except (NameError, IndexError, AttributeError) as error:
        print (error)
        print ("problems")
        flag_isTuple = False
    #
#







# ##########################
# ######  MAIN/TESTS  ######
# ##########################

def test1():
    
    r = TreeBuilder([3], [Node([2]),Node([4]),Node([5])])
    print(r.tree.children[0].name)
    
    r.tree.addToAllChildren([25,14])
    print(r.tree.children[1].name)
    print(r.tree.children[2].name)
    r.tree.children += [Node([8])]
    print(r.tree.children[3].name)
    r.tree.children[3].children = [Node([255])]
    print ("number of children in level 1: ", len(r.tree.children))
    print ("maximum height: ", r.tree.height())
    
    print ("preorder: ")
    v = r.tree.preorder()
    for i in v:
        print (i)
    #
    
    print ("disp_levels:")
    print (r.tree.disp_levels(0))
    v = NodeAnalysis()
    v.structLayout(r.tree)
    print ("a level analysis")
    print (v.disp_structLayout())
    
#


def test2(cl):
    print ("----------------------------")
    print ("now for what I want to do...")
    print (cl)
    
    a = [0]
    if len(cl) > 1:
        a = []
        for i in range(1, len(cl)):
            a += [int(cl[i])]
        #
        print (a)
    #
    t = TreeBuilder()
    t.tree.name = [3]
    t.tree.children = [Node(), Node(), Node()]
    t.tree.children[0].name = [4]
    t.tree.children[0].children = [Node()]
    t.tree.children[0].children[0].name = [5]
    t.tree.children[1].name = [40]
    t.tree.children[1].children = [Node(), Node()]
    t.tree.children[1].children[0].name = [6]
    t.tree.children[1].children[1].name = [7]
    t.tree.children[1].children[1].children = [Node()]
    t.tree.children[1].children[1].children[0].name = [8]
    t.tree.children[2].name = [50]
    
    print ("exist %s?: %s -> %s" % (a, t.existsNode(a), t.readNode(a)))
    
    print ("disp_levels:")
    print (t.tree.disp_levels(0))
    print ("maximum height: ", t.tree.height())
    
    print ("NodeAnalysis:")
    print ("a level analysis")
    v = NodeAnalysis()
    v.structLayout(t.tree)
    print (v.disp_structLayout())

    u = TreeBuilder([0])
    u.tree.children = [Node([10]), Node([11])]
    u.tree.children[1].children = [Node([210]), Node([211]), Node([212])]
    a1 = [1, 1]
    print ("try ", a1, ", depth = ", u.howDeep(a1))
    # sys.exit(0)
    a2 = [3, 1, 2, 0]
    print ("try ", a2, ", depth = ", u.howDeep(a2))
    a3 = [2, 1, 2]
    print ("try ", a3, ", depth = ", u.howDeep(a3))
    u.add_newNode(a3, [Node([3120]), Node([3121])])
    a4 = [3, 1, 2, 1]
    u.add_newNode(a4, [Node([3120]), Node([3121])])
    a5 = [4, 1, 2, 0, 0]
    u.add_newNode(a5, [Node([3120]), Node([3121])])
    u.add_newNode(a2, [Node([3120]), Node([3121])])
    w = NodeAnalysis()
    # print (u.tree)
    w.structLayout(u.tree)
    print (w.disp_structLayout())
    print ("copy current tree")
    new_u = u.copyTree()
    w.structLayout(u.tree)
    print (w.disp_structLayout())
    print ("number of nodes:      ", u.numNodes())
    print ("number of leaf nodes: ", u.numLeafNodes())
    
    print ("delete a branch ", a4)
    print ("%s exists? %s" % (a4, u.existsNode(a4)))
    
    u.deleteChild(a4)
    # print (u.tree)
    w.structLayout(u.tree)
    print (w.disp_structLayout())
    u.deleteChild(a2)
    w.structLayout(u.tree)
    print (w.disp_structLayout())
#    


# Main
if __name__ == '__main__':
    if TEST1:
        test1()
    elif TEST2:
        test2(sys.argv)
    else:
        print ("ain't got nothin to do")
    #
#
