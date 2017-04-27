#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'
'''
LinkedBinaryTree():
      __________________________________________________________
     |                                                          |      
     |    | Position  |      __Node      |   Position           |      |
     |    |           |                  |                      |      |
    \|/   |           |     | -> element |                      |      |
    |_|   |---> p --->| |_|-| -> parent  |---> position of parent      |       
          |           |     | -> left    |---> position of left child  | --->|_| 
          |           |     | -> right   |---> postiion of right child | --->|_| 
#parent                       child                                       child of child

    __Node:
        element, parent, left, right: varibles(type pointer, Position, Position, Position)
        getElement(), getParent(), getLeft(), getRight(): return value of varibles
        setElement, setParent(), setLeft(), setRight(): set value of varibles
    Position:
        container, element: varibles(type: pointer, ref)
        getContainer(), getElement(): return value of varibles

    __validate(p): Position to __Node
    __make_position(node): __Node to Position
    left(p), right(p), parent(p): Position to Position
'''
class Marker(object):
    '''To store in element of node. 
    Reduce the compute complex
    '''
    __slots__ = '__depth', '__value', '__dir', '__num', '__cross', '__clean'
    def __init__(self, depth, value, dir, num=0, cross=0, clean=0):
        self.__depth = depth
        self.__value = value
        self.__dir = dir # 0 is left, 1 is right
        self.__num = num # link num
        self.__cross = cross # cross over value
        self.__clean = clean # 0 is ok, 1 is need to be clean

    def getDepth(self):
        return self.__depth

    def getValue(self):
        return self.__value

    def getDir(self):
        return self.__dir

    def getNum(self):
        return self.__num

    def getCross(self):
        return self.__cross

    def getClean(self):
        return self.__clean

    def setValue(self, v):
        self.__value = v

    def setDepth(self, d):
        self.__depth = d

    def setNum(self, n):
        self.__num = n

    def setCross(self, c):
        self.__cross = c

    def setClean(self, l):
        self.__clean = l

class LinkedBinaryTree(object):
    '''Linked representeation of a binary tree structure.'''
    class __Node(object):
        __slots__ = '__element', '__parent', '__left', '__right'
        def __init__(self, element, parent=None, left=None, right=None):
            self.__element = element
            self.__parent = parent
            self.__left = left
            self.__right = right

        def getElement(self):
            return self.__element

        def getParent(self):
            return self.__parent

        def getLeft(self):
            return self.__left

        def getRight(self):
            return self.__right

        def setElement(self, element):
            self.__element = element

        def setParent(self, parent):
            self.__parent = parent

        def setLeft(self, left):
            self.__left = left

        def setRight(self, right):
            self.__right = right

    class Position(object):
        '''An abstrction representing the location of a single elemenet.'''
        def __init__(self, container, node):
            '''Construtor should not be invoked by user.'''
            self.__container = container # container is the tree itself. to avoid other tree's position instance
            self.__node = node

        def getElement(self):
            '''Return the element stored at this Position.'''
            return self.__node.getElement()

        def getContainer(self):
            '''Return the container of Position.'''
            return self.__container

        def getNode(self):
            return self.__node

        def __eq__(self, other):
            '''Return True if other Position represents the same location.'''
            return type(other) is type(self) and other.__node is self.__node

        def __ne__(self, other):
            '''Return True if other dose not represent the same location.'''
            return not(self == other)

    def __validate(self, p):
        '''Return associated node, if position is valid.'''
        if not isinstance(p, self.Position):
            raise TypeError('p must be proper Position type.', type(p))
        if p.getContainer() is not self:
            raise ValueError('p does not belong to this container.')
        if p.getNode().getParent() is p.getNode():
            raise ValueError('p is no longer valid.')
        return p.getNode()

    def __make_position(self, node):
        '''Return Position instance for given node (or None if no node).'''
        if node is not None: 
            return self.Position(self, node) 
        else:
            return None

    def __init__(self):
        '''Create an initially empty binary tree.'''
        self.__root = None
        self.__size = 0

    def root(self):
        '''Return Position resenting the tree's root(or None if empty).'''
        return self.__make_position(self.__root)

    def __len__(self):
        '''Return the total number of elements in the tree.'''
        return self.__size

    def parent(self, p):
        '''Returen Position representing p's parent (or None if p is root).'''
        node = self.__validate(p)
        return self.__make_position(node.getParent())

    def left(self, p):
        '''Return the Position of p's left child (or None if no left child).'''
        node = self.__validate(p)
        return self.__make_position(node.getLeft())

    def right(self, p):
        '''Return the Position of p's right child (or None if no left child).'''
        node = self.__validate(p)
        return self.__make_position(node.getRight())

    def num_children(self, p):
        '''Return the number of children that Position P has.'''
        node = self.__validate(p)
        count = 0
        if node.getLeft() is not None:
            count += 1
        if node.getRight() is not None:
            count += 1
        return count

    def add_root(self, e):
        '''Place element e at the root of an empty tree and return noe Position.
        Raise ValueError if tree nonempty
        '''
        if self.__root is not None:
            raise ValueError('Root exists.')
        self.__size = 1
        self.__root = self.__Node(e)
        return self.__make_position(self.__root)

    def add_left(self, p, e):
        '''creat a new left child for position p, storing element e.
        return the position of new node.
        raise valueerror if position p is invalid or p already has a left child.
        '''
        node = self.__validate(p)
        if node.getLeft() is not None:
            raise valueerror('left child exists.')
        self.__size += 1
        node.setLeft(self.__Node(e, node))
        return self.__make_position(node.getLeft())

    def add_right(self, p, e):
        '''Creat a new right child for Position p, storing element e.
        Return the Position of new node.
        Raise ValueError if Position p is invalid or p already has a right child.
        '''
        node = self.__validate(p)
        if node.getRight() is not None:
            raise ValueError('Left child exists.')
        self.__size += 1
        node.setRight(self.__Node(e, node))
        return self.__make_position(node.getRight())

    def replace(self, p, e):
        '''Replace the element at position p with e, and return old element.'''
        node = self.__validate(p)
        old = node.getElement()
        node.setElement(e)
        return old

    def attach(self, p, t1, t2):
        '''Attach tree t1 an t2 as left and right subtrees of external p.'''
        node = self.__validate(p)
        if not self.is_leaf(p):
            raise ValueError('position must be leaf.')
        if not type(self) is type(t1) is type(t2):
            raise TypeError('Tree types must match.')
        self.__size += len(t1) + len(t2)
        if not t1.is_empty():
            t1.__root.setParent(node)
            node.setLeft(t1.__root)
            t1.__root = None
            t1.__size = 0
        if not t2.is_empty():
            t2.__root.setParent(node)
            node.setRight(t2.__root)
            t2.__root = None
            t2.__size = 0

    def is_root(self, p):
        '''Return True if Position p represents the root of the tree.'''
        return self.root() == p

    def is_leaf(self, p):
        '''Return True if Position p does not have any children.'''
        return self.num_children(p) == 0

    def is_empty(self):
        '''Return True if the tree is empty.'''
        return len(self) == 0

    def sibling(self, p):
        '''Retrun a Position representing p's sibling (or None if no sibling).'''
        parent = self.parent(p)
        if parent is None:
            return None
        else:
            if p == self.left(parent):
                return self.right(parent)
            else:
                return self.left(parent)
    
    def children(self, p):
        '''Generate an iteration of Positions represnting p's children.'''
        if self.left(p) is not None:
            yield self.left(p)
        if self.right(p) is not None:
            yield self.right(p)

    def depth(self, p):
        '''Return the number of levels separation Position p from the root.'''
        if self.is_root(p):
            return 0
        else:
            return 1 + self.depth(self.parent(p))

    def height(self, p=None):
        '''Return the height of the subtree rooted at Position p.
        If p is None return the height of entire tree.
        '''
        if p is None:
            p = self.root()
        if self.is_leaf(p):
            return 0
        return 1 + max(self. height(c) for c in self.children(p))

    def delete(self, p):
        '''Delete the node at Position p, and replace it with its child, if any.
        Return the element that had been storedat Postion p.
        Raise ValueError if Position p is invalid or p has two children.
        '''
        node = self.__validate(p)
        if self.num_children(p) == 2:
            raise ValueError('p has two children tree.')
        child = node.getLeft() if node.getLeft() is not None else node.getRight()
        if child is not None:
            child.setParent(node.parent)
        if node is self.__root:
            self.__root = child
        else:
            parent = node.getParent()
            if node is parent.getLeft():
                parent.setLeft(child)
            else:
                parent.setRight(child)
        self.__size -= 1
        node.setParent(node)
        return node.getElement()

    def __iter__(self):
        '''Generate an iteration of tree's elements.'''
        for p in self.inorder():
            yield p.getElement()

    def inorder(self):
        '''Generate a preorder iteration of positions in the tree.'''
        if not self.is_empty():
            for p in self.__subtree_inorder(self.root()):
                yield p

    def __subtree_preorder(self,p):
        '''Generate a preorder iteration of positions in subtree rooted at p.'''
        yield p
        for c in self.children(p):
            for other in self.__subtree_preorder(c):
                yield other

    def __subtree_postorder(self, p):
        '''Generate a postorder iteration of positions in subtree rooted at P.'''
        for c in self.children(p):
            for other in self.__subtree_postorder(c):
                yield other
        yield p

    def __subtree_inorder(self, p):
        '''Generate a inorder iteration of positions in subtree rooted at p.(only for binary tree.'''
        if self.left(p) is not None:
            for other in self.__subtree_inorder(self.left(p)):
                yield other
        yield p
        if self.right(p) is not None:
            for other in self.__subtree_inorder(self.right(p)):
                yield other

    def inorder_indent(self, p, d=0):
        '''Print inorder representation of subtree of T rooted at p at start with 2d*blank.
        Develper's tools
        '''
        if not self.is_empty():
            for other in self.__subtree_inorder(p):
                node = self.__validate(other)
                dep = node.getElement().getDepth()  # depth of node
                dir = node.getElement().getDir() # direction of node 
                num = node.getElement().getNum() # direction of node 
                print(2*dep*' '+ str(dir) + ' ' + str(num) )

    def preorder_indent(self, p, d=0):
        node = self.__validate(p)
        dep = d + node.getElement().getDepth()  # depth of node
        dir = node.getElement().getDir() # direction of node 
        num = node.getElement().getNum() # direction of node 
        print(2*dep*' '+ str(dir) + ' -> ' + str(num) )
        for c in self.children(p):
            self.preorder_indent(c, d)

    def delete_subtree(self, p):
        '''Delete the node at Position p, and its child, if any.
        Only work in post order.
        '''
        if not self.is_empty():
            for other in self.__subtree_postorder(p):
                node = self.__validate(other)
                if not self.is_root(p):
                    if other == p and not self.is_root(p): # root of subtree
                        p_p = self.parent(other) # position of parent
                        p_n = self.__validate(p_p) # node of parent
                        if self.left(p_p) == other:
                            p_n.setLeft(None)
                        elif self.right(p_p) == other:
                            p_n.setRight(None)
                        else:
                            raise ValueError('Neighter of left or right child')
                    else:
                        node.setParent(None)
                node.setElement(None)
                node.setLeft(None)
                node.setRight(None)
                other = None

    def pruning(self, p, heter_snp, level_pos, tree_p):
        '''Decide the SNP linkage result(based on element value).
        And crop non SNP linkage branch and subtree.
        And return the root of tree or the break point(new root of tree)
        '''
        pruning_f = open(tree_p, 'a') 
        if not self.is_empty():
            for other in self.__subtree_preorder(p):
                node = self.__validate(other)
                dep = node.getElement().getDepth()  # left child element value
                if dep == 0: # root
                    continue
                if self.num_children(other) == 0: # leaf
                    continue
                left_p  = self.left(other)  # left child position
                right_p = self.right(other) # right child position
                if left_p is not None and right_p is not None:
                    left_num = self.__validate(left_p).getElement().getNum()  # left child element link num
                    right_num = self.__validate(right_p).getElement().getNum()  # right child element link num
                    left_value = self.__validate(left_p).getElement().getValue()  # left child element value
                    right_value = self.__validate(right_p).getElement().getValue()  # right child element value
                    left_cross = self.__validate(left_p).getElement().getCross() # value + crossover
                    right_cross = self.__validate(right_p).getElement().getCross() # value + crossover

                    ind = level_pos[dep+1]
                    heter_base = heter_snp[ind] # heter_snp index of node = dep-1; heter_base: heter snp bases
                    #print('%d\t%d' % (left_num, right_num))
                    num_dir, cross_dir = 0, 0 # link num direction result, crossover link num direction result 
                    if right_num > left_num:
                        num_dir = 1
                    if right_cross > left_cross:
                        cross_dir = 1
                    com = (num_dir == cross_dir) # compare link num and crossover 
                    sig_num = self.significant(left_num, right_num) # test two child if significant
                    sig_cross = self.significant(left_cross, right_cross) # test two child if significant
                    AE = self.align_error(num_dir, cross_dir, sig_num, sig_cross, heter_base) # if it is alignment error

                    #print(dep, left_num, right_num, left_cross, right_cross, AE)
                    if com is True and AE is not True:
                        pruning_f.write('%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (int(100*left_value),int(100*right_value),left_num,right_num,left_cross,right_cross, max(right_num, left_num), min(right_num, left_num)))
                    #sig = True
                    if left_num is not None and right_num is not None:
                        if left_num > right_num: # left child element value is larger, delete right child tree
                            self.delete_subtree(right_p)
                            if com is not True:
                                self.__validate(left_p).getElement().setClean(2)  # left child element need to clean as homo snp
                            if AE is True:
                                self.__validate(left_p).getElement().setClean(1)  # left child element need to clean as seq/align error
                        elif left_num < right_num: # right child element value is larger, delete left child tree
                            self.delete_subtree(left_p)
                            if com is not True:
                                self.__validate(right_p).getElement().setClean(2)  # right child element need to clean as homo snp
                            if AE is True:
                                self.__validate(right_p).getElement().setClean(1)  # right child element need to clean as seq/align error
                        elif left_num == right_num == 0: # no link information 
                            if left_cross > right_cross: # use cross information
                                self.delete_subtree(right_p)
                            elif right_cross > left_cross:
                                self.delete_subtree(left_p)
                            else: # no link info and crossover info
                                raise ValueError('depth', dep, 'Should be break point.')
                        else: # right == left and != 0
                            print('Wrong at', dep,'-', dep+1, left_num, '=', right_num)
                            return dep + 1

    def align_error(self, num_dir, cross_dir, sig_num, sig_cross, heter_base):
        node_base = heter_base[num_dir]
        cross_base = heter_base[cross_dir]
        if cross_base == 'R' and sig_cross is True and sig_num is not True:
            return True
        if node_base == 'R' and sig_num is True and sig_cross is not True:
            return True
        return False

    def significant(self, v1, v2):
        '''Return True if SNP(same depth) is significant'''
        min_num = min(v1, v2)
        if min_num == 0:
            return True
        return abs(v1-v2)/min_num >= 0.55

    def add_value_left(self, d, v, direct):
        '''Add value=v to all node element in depth=d.
           dir == 1 is right
           dir = 0 is left 
        '''
        if not self.is_empty():
            for other in self.__subtree_preorder(self.root()):
                node = self.__validate(other)
                dep = node.getElement().getDepth()  # depth of node
                dir = node.getElement().getDir() # direction of node 
                if dep == d and dir == direct: # make sure parent's  is same as dir
                    left_c = self.left(other)
                    if left_c is not None:
                        node = self.__validate(left_c)
                        mar = node.getElement()
                        prev_v = mar.getValue()
                        mar.setValue(mar.getValue() + v)
                        mar.setNum(mar.getNum() + 1)
                        assert prev_v != mar.getValue(), 'Value add error'

    def add_value_right(self, d, v, direct):
        '''Add value=v to all node element in depth=d.
           dir == 1 is right
           dir = 0 is left 
        '''
        if not self.is_empty():
            for other in self.__subtree_postorder(self.root()):
                node = self.__validate(other)
                dep = node.getElement().getDepth()  # depth of node
                dir = node.getElement().getDir() # direction of node 
                if dep == d and dir == direct: # make sure parent's  is same as dir
                    right_c = self.right(other) # left child position
                    if right_c is not None:
                        node = self.__validate(right_c)
                        mar = node.getElement()
                        prev_v = mar.getValue()
                        mar.setValue(mar.getValue() + v)
                        mar.setNum(mar.getNum() + 1)
                        assert prev_v != node.getElement().getValue(), 'Value add error'

    def add_cross_left(self, d, c, direct):
        '''Add cross=v to all node element in depth=d.
        '''
        if not self.is_empty():
            for other in self.__subtree_preorder(self.root()):
                node = self.__validate(other)
                dep = node.getElement().getDepth()  # depth of node
                dir = node.getElement().getDir() # direction of node 
                if dep == d and dir == direct: # make sure parent's  is same as dir
                    left_c = self.left(other)
                    if left_c is not None:
                        mar = self.__validate(left_c).getElement()
                        prev_c = mar.getCross()
                        mar.setCross(prev_c + c)
                        assert prev_c != mar.getCross(), 'Value add error'

    def add_cross_right(self, d, c, direct):
        '''Add cross=v to all node element in depth=d.
        '''
        if not self.is_empty():
            for other in self.__subtree_preorder(self.root()):
                node = self.__validate(other)
                dep = node.getElement().getDepth()  # depth of node
                dir = node.getElement().getDir() # direction of node 
                if dep == d and dir == direct: # make sure parent's  is same as dir
                    right_c = self.right(other)
                    if right_c is not None:
                        mar = self.__validate(right_c).getElement()
                        prev_c = mar.getCross()
                        mar.setCross(prev_c + c)
                        assert prev_c != mar.getCross(), 'Value add error'

    def linkage_result(self):
        '''Get conclusion of heter-snp-marker linkage infomation.'''
        h = self.height()
        sub_l = []*h #left
        sub_r = []*h #right
        for other in self.__subtree_preorder(self.left(self.root())):# sub left tree
            node = self.__validate(other)
            dep = node.getElement().getDepth()  # depth of node
            dir = node.getElement().getDir() # direction of node 
            sub_l.append(dir)
        for other in self.__subtree_preorder(self.right(self.root())):# sub right tree
            node = self.__validate(other)
            dep = node.getElement().getDepth()  # depth of node
            dir = node.getElement().getDir() # direction of node 
            sub_r.append(dir)
        return sub_l, sub_r

    def setdefault(self, d, v):
        '''Init 2 child of depth d and set default v.
        Traverse tree through postorder.
        '''
        v = float(v)
        before_v = 0.0 # set 0 before depth
        for other in self.__subtree_preorder(self.root()):
            #dep = self.__validate(other).getElement().getDepth()  # depth of node
            dep = self.depth(other)
            if dep < d-1 and self.num_children(other) == 0:
                self.add_left(other,  Marker(dep+1, before_v, 0))
                self.add_right(other, Marker(dep+1, before_v, 1))
            elif dep == d-1 and self.num_children(other) == 0:
                self.add_left(other,  Marker(dep+1, v, 0))
                self.add_right(other, Marker(dep+1, v, 1))

    def clean(self, level_s):
        '''Clean the tree, rm seq/align error and homo snp'''
        mis_level = None
        clean_type = 0
        for other in self.__subtree_preorder(self.root()):
            mar = self.__validate(other).getElement()
            dep = mar.getDepth()  # depth of mar
            clean_type = mar.getClean() # clean value: 0=do not clean 1=seq/align error 2=homo snp
            mis_level = dep 
            if clean_type == 1:
                self.delete_depth(dep)
                break
            if clean_type == 2:
                self.delete_depth(dep)
                break
        if clean_type == 1:
            print('    rm seq-error/mis-aligned SNP at level:%d' % dep)
        elif clean_type == 2:
            print('    rm homo/ambiguous-heter SNP at level:%d' % dep)
        else:
            mis_level = None
        return mis_level

    def delete_depth(self, d):
        for other in self.__subtree_postorder(self.root()):
            mar = self.__validate(other).getElement()
            dep = mar.getDepth()  # depth of mar
            if dep >= d:
                self.delete_subtree(other)

if __name__ == '__main__':
    import cpuCheck, os
    pid = os.getpid()
    print('origin:', cpuCheck.mem(pid))
    
    # Test for binary tree init and print out
    lbt = LinkedBinaryTree()
    p0 = lbt.add_root('/')
    p00 = lbt.add_left(p0,'Document')
    pt = lbt.left(p0)
    print('test: left():', pt ==p00)
    p000 = lbt.add_left(p00,'PDFs')
    p001 = lbt.add_right(p00,'PNGs')
    p01 = lbt.add_right(p0,'Desktop')
    p010 = lbt.add_left(p01,'Wallpapers')
    p011 = lbt.add_right(p01,'Folds')
    #for x in lbt:
    #    print(x)
    lbt.inorder_indent(lbt.root(), 0)

    # Test for memory realse
    lbt = LinkedBinaryTree()
    p = lbt.add_root('root')
    p0 = lbt.add_left(p,'left')
    p1 = lbt.add_right(p,'right')
    for x in range(0,100):
        p0 = lbt.add_left(p0, [0]*1000)
    for x in range(0,100):
        p1 = lbt.add_right(p1, [1]*1000)
    print('init tree:', cpuCheck.mem(pid))
    p0 = lbt.left(p)
    p0 = lbt.left(p0)
    p1 = lbt.right(p)
    p1 = lbt.right(p1)
    lbt.delete_subtree(p0)
    lbt.delete_subtree(p1)
    lbt.inorder_indent(lbt.root())
    print('delete tree:', cpuCheck.mem(pid))
    
    # Test for add_value, pruning, break and linkage result
    lbt = LinkedBinaryTree()
    p0 = lbt.add_root('root')
    lbt.setdefault(1,1)
    lbt.setdefault(5,0)
    lbt.inorder_indent(lbt.root())
    # 0 1 0 1 1 *3
    lbt.add_value_right(1, 3, 0)
    lbt.add_value_right(1, 2.5, 1) # mistakejZ
    lbt.add_value_left(2, 3, 1)
    lbt.add_value_right(3, 3, 0)
    lbt.add_value_right(3, 2.5, 1) # mistake
    lbt.add_value_right(4, 3, 1)
    # 1 0 1 0 0 *3
    '''
    left tree
    0 1 0 1 1 
    right tree
    1 0 1 0 0
    '''
    lbt.add_value_left(1, 3, 1)
    lbt.add_value_right(2, 3, 0)
    lbt.add_value_right(2, 2, 1) # mistake ai
    lbt.add_value_left(3, 3, 1)
    lbt.add_value_left(4, 3, 0)
    lbt.pruning(lbt.root(),5) # depth 3
    #lbt.inorder_indent(lbt.root())
    left_t = lbt.left(lbt.root())
    right_t = lbt.right(lbt.root())
    lbt.preorder_indent(left_t)
    lbt.preorder_indent(right_t)
    lbt.linkage_result()

    '''
    lbt.inorder_indent(lbt.root(), 0)
    break_point = lbt.pruning(lbt.root(),0) # depth 3
    print('pruning leve 0')
    lbt.inorder_indent(lbt.root())
    break_point = lbt.pruning(lbt.root(),1) # depth 3
    print('pruning leve 1')
    lbt.inorder_indent(lbt.root())
    break_point = lbt.pruning(lbt.root(),2) # depth 3
    print('pruning leve 2')
    lbt.inorder_indent(lbt.root())
    break_point = lbt.pruning(lbt.root(),3) # depth 3
    print('pruning leve 3')
    lbt.inorder_indent(lbt.root())
    '''
