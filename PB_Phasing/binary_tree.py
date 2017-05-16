#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
LinkedBinaryTree():
      __________________________________________________________
     |                                                          |
     |    | Position  |      _Node      |   Position           |      |
     |    |           |                  |                      |      |
    \|/   |           |     | -> element |                      |      |
    |_|   |---> p --->| |_|-| -> parent  |---> position of parent      |
          |           |     | -> left    |---> position of left child  | --->|_|
          |           |     | -> right   |---> postiion of right child | --->|_|
#parent                       child                                       child of child

    _Node:
        element, parent, left, right: varibles(type pointer, Position, Position, Position)
        get_element(), get_parent(), get_left(), get_right(): return value of varibles
        set_element, set_parent(), set_left(), set_right(): set value of varibles
    Position:
        container, element: varibles(type: pointer, ref)
        get_container(), get_element(): return value of varibles

    __validate(p): Position to _Node
    __make_position(node): _Node to Position
    left(p), right(p), parent(p): Position to Position
'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

class Marker(object):
    '''To store in element of node. Reduce the compute complex.'''

    __slots__ = '__depth', '__num', '__direct', '__cross', '__clean', '__significant'
    def __init__(self, depth, num, direct, cross=0, clean=0):
        self.__depth = depth
        self.__direct = direct # 0 is left, 1 is right
        self.__num = num # link num
        self.__cross = cross # cross over value
        self.__clean = clean # 0 is ok, 1 is need to be clean
        self.__significant = True

    def get_depth(self):
        '''Return depth information of the class Marker'''
        return self.__depth

    def get_dir(self):
        '''Retrun dirction information of the class Marker'''
        return self.__direct

    def get_num(self):
        '''Return link number information of the class Marker'''
        return self.__num

    def get_cross(self):
        '''Return crossover informatino info of the class Marker'''
        return self.__cross

    def get_clean(self):
        '''Return if need to be clean informatino info of the class Marker'''
        return self.__clean

    def get_sig(self):
        '''Return if linkage number is significant.'''
        return self.__significant

    def set_depth(self, d):
        '''Set depth information of the class Marker'''
        self.__depth = d

    def set_num(self, n):
        '''Set link number information of the class Marker'''
        self.__num = n

    def set_cross(self, c):
        '''Set crossover informatino info of the class Marker'''
        self.__cross = c

    def set_clean(self, l):
        '''Set if need to be clean informatino info of the class Marker'''
        self.__clean = l

    def set_sig(self, s):
        '''Set if linkage number is sinificant.'''
        self.__significant = s

    def delete(self):
        '''Set all term as None.'''
        self.__depth = None
        self.__direct = None
        self.__num = None
        self.__cross = None
        self.__clean = None
        self.__significant = None

class LinkedBinaryTree(object):
    '''Linked representeation of a binary tree structure.'''
    class _Node(object):
        __slots__ = '__element', '__parent', '__left', '__right'
        def __init__(self, element, parent=None, left=None, right=None):
            self.__element = element
            self.__parent = parent
            self.__left = left
            self.__right = right

        def get_element(self):
            '''Return element.'''
            return self.__element

        def get_parent(self):
            '''Return parent Position.'''
            return self.__parent

        def get_left(self):
            '''Return left child Position.'''
            return self.__left

        def get_right(self):
            '''Return right child Position.'''
            return self.__right

        def set_element(self, element):
            '''Set element.'''
            self.__element = element

        def set_parent(self, parent):
            '''Set Parent.'''
            self.__parent = parent

        def set_left(self, left):
            '''Set left child.'''
            self.__left = left

        def set_right(self, right):
            '''Set right child.'''
            self.__right = right

    class Position(object):
        '''An abstrction representing the location of a single elemenet.'''
        def __init__(self, container, node):
            '''Construtor should not be invoked by user.'''
            self.__container = container # container is the tree itself.
                                         # To avoid other tree's position instance
            self.__node = node

        def get_element(self):
            '''Return the element stored at this Position.'''
            return self.__node.get_element()

        def get_container(self):
            '''Return the container of Position.'''
            return self.__container

        def get_node(self):
            '''Return the node (element) of Position.'''
            return self.__node

        def __eq__(self, other):
            '''Return True if other Position represents the same location.'''
            return type(other) is type(self) and other.get_node() is self.get_node()

        def __ne__(self, other):
            '''Return True if other dose not represent the same location.'''
            return not self == other

    def __validate(self, p):
        '''Return associated node, if position is valid.'''
        if not isinstance(p, self.Position):
            raise TypeError('p must be proper Position type.', type(p))
        if p.get_container() is not self:
            raise ValueError('p does not belong to this container.')
        if p.get_node().get_parent() is p.get_node():
            raise ValueError('p is no longer valid.')
        return p.get_node()

    def __make_position(self, node):
        '''Return Position instance for given node (or None if no node).'''
        if node is not None:
            return self.Position(self, node)
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

    def set_size(self, n=0):
        '''Set the total number of elements in the tree.'''
        self.__size = n

    def parent(self, p):
        '''Returen Position representing p's parent (or None if p is root).'''
        node = self.__validate(p)
        return self.__make_position(node.get_parent())

    def left(self, p):
        '''Return the Position of p's left child (or None if no left child).'''
        node = self.__validate(p)
        return self.__make_position(node.get_left())

    def right(self, p):
        '''Return the Position of p's right child (or None if no left child).'''
        node = self.__validate(p)
        return self.__make_position(node.get_right())

    def num_children(self, p):
        '''Return the number of children that Position P has.'''
        node = self.__validate(p)
        count = 0
        if node.get_left() is not None:
            count += 1
        if node.get_right() is not None:
            count += 1
        return count

    def add_root(self, e):
        '''Place element e at the root of an empty tree and return none Position.
        Raise ValueError if tree nonempty.'''
        if self.__root is not None:
            raise ValueError('Root exists.')
        self.__size = 1
        self.__root = self._Node(e)
        return self.__make_position(self.__root)

    def set_root(self, e):
        '''Set root of a none empty tree'''
        if self.__root is not None:
            self.__root = self._Node(e)
            return self.__make_position(self.__root)
        return None

    def add_left(self, p, e):
        '''creat a new left child for position p, storing element e.
        return the position of new node.
        raise valueerror if position p is invalid or p already has a left child.
        '''
        node = self.__validate(p)
        if node.get_left() is not None:
            raise ValueError('left child exists.')
        self.__size += 1
        node.set_left(self._Node(e, node))
        return self.__make_position(node.get_left())

    def add_right(self, p, e):
        '''Creat a new right child for Position p, storing element e.
        Return the Position of new node.
        Raise ValueError if Position p is invalid or p already has a right child.
        '''
        node = self.__validate(p)
        if node.get_right() is not None:
            raise ValueError('Left child exists.')
        self.__size += 1
        node.set_right(self._Node(e, node))
        return self.__make_position(node.get_right())

    def replace(self, p, e):
        '''Replace the element at position p with e, and return old element.'''
        node = self.__validate(p)
        old = node.get_element()
        node.set_element(e)
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
            t1.root().set_parent(node)
            node.set_left(t1.root())
            t1.set_root(None)
            t1.set_size(0)
        if not t2.is_empty():
            t2.root().set_parent(node)
            node.set_right(t2.root())
            t2.set_root(None)
            t2.set_size(0)

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
        if p == self.left(parent):
            return self.right(parent)
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
        child = node.get_left() if node.get_left() is not None else node.get_right()
        if child is not None:
            child.set_parent(node.parent)
        if node is self.__root:
            self.__root = child
        else:
            parent = node.get_parent()
            if node is parent.get_left():
                parent.set_left(child)
            else:
                parent.set_right(child)
        self.__size -= 1
        node.set_parent(node)
        return node.get_element()

    def __iter__(self):
        '''Generate an iteration of tree's elements.'''
        for p in self.inorder():
            yield p.get_element()

    def inorder(self):
        '''Generate a preorder iteration of positions in the tree.'''
        if not self.is_empty():
            for p in self.__subtree_inorder(self.root()):
                yield p

    def __subtree_preorder(self, p):
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
        ''' Generate a inorder iteration of positions in subtree rooted at p.
            (only for binary tree)'''
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
                dep = node.get_element().get_depth()  # depth of node
                direct = node.get_element().get_dir() # directection of node
                num = node.get_element().get_num() # directection of node
                print(2*dep*' '+ str(direct) + ' ' + str(num))

    def preorder_indent(self, p, d=0):
        node = self.__validate(p)
        dep = d + node.get_element().get_depth()  # depth of node
        direct = node.get_element().get_dir() # directection of node
        num = node.get_element().get_num() # directection of node
        print(2*dep*' '+ str(direct) + ' -> ' + str(num))
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
                            p_n.set_left(None)
                        elif self.right(p_p) == other:
                            p_n.set_right(None)
                        else:
                            raise ValueError('Neighter of left or right child')
                    else:
                        node.set_parent(None)
                node.get_element().delete() # delete item in class Marker
                node.set_element(None)
                node.set_left(None)
                node.set_right(None)
                other = None

    def pruning(self, p, tree_p):
        '''Decide the SNP linkage result(based on element value).
        And crop non SNP linkage branch and subtree.
        And return the root of tree or the break point(new root of tree)
        '''
        pruning_f = open(tree_p, 'a')
        result = [] # result list: tree node Position
        if not self.is_empty():
            for other in self.__subtree_preorder(p):
                node = self.__validate(other)
                #print(other, node.get_element())
                dep = node.get_element().get_depth()  # left child element value
                if dep == 0: # root
                    continue
                if self.num_children(other) == 0: # leaf
                    continue
                left_p = self.left(other)  # left child position
                right_p = self.right(other) # right child position
                if left_p is not None and right_p is not None:
                    # left and right child element (class Marker)
                    left_mar = self.__validate(left_p).get_element()
                    right_mar = self.__validate(right_p).get_element()
                    left_num = left_mar.get_num()  # left child element link num
                    right_num = right_mar.get_num()  # right child element link num
                    left_cross = left_mar.get_cross() # value + crossover
                    right_cross = right_mar.get_cross() # value + crossover

                    '''homo snp removing.'''
                    # link num directection result, crossover link num directection result
                    num_direct, cross_direct = 0, 0
                    if right_num > left_num:
                        num_direct = 1
                    if right_cross > left_cross:
                        cross_direct = 1
                    com = (num_direct == cross_direct) # compare link num and crossover direction

                    '''linkage number and crossover significant check.'''
                    # test two children's linkage number if significant
                    sig_num = self.significant(left_num, right_num)
                    # test two children's crossover if significant
                    sig_cross = self.significant(left_cross, right_cross)

                    # if it is alignment error
                    ae = sig_num is not True and sig_cross is not True
                    if com is True and ae is not True:
                        pruning_f.write('%d\t%d\t%d\t%d\t%d\t%d\n' % \
                                        (left_num, right_num, left_cross, right_cross, \
                                         max(right_num, left_num), min(right_num, left_num)))
                    if left_num is not None and right_num is not None:
                        # left child element value is larger, delete right child tree
                        if left_num >= right_num:
                            result.append(left_p) # add left child node Position to result list
                            self.delete_subtree(right_p)
                            if com is not True:
                                # left child element need to clean as homo snp
                                self.__validate(left_p).get_element().set_clean(2)
                            if ae is True:
                                # left child element need to clean as seq/align error
                                self.__validate(left_p).get_element().set_clean(1)
                        # right child element value is larger, delete left child tree
                        elif left_num < right_num:
                            result.append(right_p) # add right child node Position to result list
                            self.delete_subtree(left_p)
                            if com is not True:
                                # right child element need to clean as homo snp
                                self.__validate(right_p).get_element().set_clean(2)
                            if ae is True:
                                # right child element need to clean as seq/align error
                                self.__validate(right_p).get_element().set_clean(1)
                        elif left_num == right_num == 0: # no link information
                            if left_cross > right_cross: # use cross information
                                result.append(left_p) # add left child node Position to result list
                                self.delete_subtree(right_p)
                            elif right_cross > left_cross:
                                result.append(right_p) # add right child node Position to result list
                                self.delete_subtree(left_p)
                            else: # no link info and crossover info
                                #raise ValueError('depth', dep, 'Should be break point.')
                                # break point
                                result.append(left_p) # add left child node Position to result list
                                self.delete_subtree(right_p) # rm right subtree (all 2 node in same level)
        return result # result list: tree node Position

    def significant(self, v1, v2):
        '''Return True if SNP(same depth) is significant'''
        min_num = min(v1, v2)
        max_num = max(v1, v2)
        if min_num == 0:
            return True
        return max_num > 2*min_num

    def add_value_left(self, p, d, n, directect):
        '''Add value=v to all node element in depth=d.
           direct == 1 is right
           direct = 0 is left
        '''
        if not self.is_empty():
            for other in self.__subtree_preorder(p):
                node = self.__validate(other)
                dep = node.get_element().get_depth()  # depth of node
                direct = node.get_element().get_dir() # directection of node
                if dep == d and direct == directect: # make sure parent's  is same as direct
                    left_c = self.left(other)
                    if left_c is not None:
                        node = self.__validate(left_c)
                        mar = node.get_element()
                        prev_n = mar.get_num()
                        mar.set_num(prev_n + n)
                        assert prev_n != mar.get_num(), 'Value add error'

    def add_value_right(self, p, d, n, directect):
        '''Add value=v to all node element in depth=d.
           direct == 1 is right
           direct = 0 is left
        '''
        if not self.is_empty():
            for other in self.__subtree_preorder(p):
                node = self.__validate(other)
                dep = node.get_element().get_depth()  # depth of node
                direct = node.get_element().get_dir() # directection of node
                if dep == d and direct == directect: # make sure parent's  is same as direct
                    right_c = self.right(other) # left child position
                    if right_c is not None:
                        node = self.__validate(right_c)
                        mar = node.get_element()
                        prev_n = mar.get_num()
                        mar.set_num(prev_n + n)
                        assert prev_n != mar.get_num(), 'Value add error'

    def add_cross_left(self, p, d, c, directect):
        '''Add cross=v to all node element in depth=d.
        '''
        if not self.is_empty():
            for other in self.__subtree_preorder(p):
                node = self.__validate(other)
                dep = node.get_element().get_depth()  # depth of node
                direct = node.get_element().get_dir() # directection of node
                if dep == d and direct == directect: # make sure parent's  is same as direct
                    left_c = self.left(other)
                    if left_c is not None:
                        mar = self.__validate(left_c).get_element()
                        prev_c = mar.get_cross()
                        mar.set_cross(prev_c + c)
                        assert prev_c != mar.get_cross(), 'Value add error'

    def add_cross_right(self, p, d, c, directect):
        '''Add cross=v to all node element in depth=d.
        '''
        if not self.is_empty():
            for other in self.__subtree_preorder(p):
                node = self.__validate(other)
                dep = node.get_element().get_depth()  # depth of node
                direct = node.get_element().get_dir() # directection of node
                if dep == d and direct == directect: # make sure parent's  is same as direct
                    right_c = self.right(other)
                    if right_c is not None:
                        mar = self.__validate(right_c).get_element()
                        prev_c = mar.get_cross()
                        mar.set_cross(prev_c + c)
                        assert prev_c != mar.get_cross(), 'Value add error'

    def linkage_result(self, l):
        '''Get conclusion of heter-snp-marker linkage infomation.
        l is length of heter_snp.
        Result list: [0/1 in heter_snp, if break point, if significant]
        '''
        sub_l = [[None, None, None] for x in range(0, l)] # list for sub left tree
        for other in self.__subtree_preorder(self.left(self.root())):# sub left tree
            node = self.__validate(other) # node in tree
            mar = node.get_element() # class Marker
            dep = mar.get_depth()  # depth of Marker
            direct = mar.get_dir() # directection of Marker
            break_point =  mar.get_num() == -1 # linkage number of Marker == -1 is break point
            significant = mar.get_sig() # significant of Marker's linkage number
            sub_l[dep-1] = [direct, break_point, significant]

        sub_r = [[None, None, None] for x in range(0, l)] # list for sub right tree
        for other in self.__subtree_preorder(self.right(self.root())):# sub right tree
            node = self.__validate(other) # node in tree
            mar = node.get_element() # class Marker
            dep = mar.get_depth()  # depth of Marker
            direct = mar.get_dir() # directection of Marker
            break_point =  mar.get_num() == -1 # linkage number of Marker == -1 is break point
            significant = mar.get_sig() # significant of Marker's linkage number
            sub_r[dep-1] = [direct, break_point, significant]
        return sub_l, sub_r

    def setdefault(self, p, d, n):
        '''Init 2 child of depth d and set default v.
        Traverse tree through postorder.
        '''
        for other in self.__subtree_preorder(p):
            dep = 0 # depth of node
            if other != self.root():
                dep = self.__validate(other).get_element().get_depth()
            if dep < d-1 and self.num_children(other) == 0:
                self.add_left(other, Marker(dep+1, 0, 0))
                self.add_right(other, Marker(dep+1, 0, 1))
            elif dep == d-1 and self.num_children(other) == 0:
                self.add_left(other, Marker(dep+1, n, 0))
                self.add_right(other, Marker(dep+1, n, 1))

    def clean(self, left_right_p):
        '''Return the min level of node need to be clean.'''
        for p in left_right_p:
            for other in self.__subtree_preorder(p):
                mar = self.__validate(other).get_element()
                dep = mar.get_depth()  # depth of mar
                clean_type = mar.get_clean() # clean value: 0=do not clean 1=seq/align error 2=homo snp
                if clean_type == 1:
                    print('    rm seq-error/mis-aligned SNP at level:%d' % dep, end=' ')
                    return dep
                elif clean_type == 2:
                    print('    rm homo/ambiguous-heter SNP at level:%d' % dep, end=' ')
                    return dep

    def delete_depth(self, p, d):
        '''Delete all tree after depth d.'''
        for other in self.__subtree_postorder(p):
            mar = self.__validate(other).get_element()
            dep = mar.get_depth()  # depth of mar
            if dep >= d:
                self.delete_subtree(other)

if __name__ == '__main__':
    import cpuCheck, os
    pid = os.getpid()
    print('origin:', cpuCheck.mem(pid))

    # Test for binary tree init and print out
    lbt = LinkedBinaryTree()
    p0 = lbt.add_root('/')
    p00 = lbt.add_left(p0, 'Document')
    pt = lbt.left(p0)
    print('test: left():', pt == p00)
    p000 = lbt.add_left(p00, 'PDFs')
    p001 = lbt.add_right(p00, 'PNGs')
    p01 = lbt.add_right(p0, 'Desktop')
    p010 = lbt.add_left(p01, 'Wallpapers')
    p011 = lbt.add_right(p01, 'Folds')
    #for x in lbt:
    #    print(x)
    lbt.inorder_indent(lbt.root(), 0)

    # Test for memory realse
    lbt = LinkedBinaryTree()
    p = lbt.add_root('root')
    p0 = lbt.add_left(p, 'left')
    p1 = lbt.add_right(p, 'right')
    for x in range(0, 100):
        p0 = lbt.add_left(p0, [0]*1000)
    for x in range(0, 100):
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
    lbt.setdefault(1, 1)
    lbt.setdefault(5, 0)
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
    lbt.pruning(lbt.root(), 5) # depth 3
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
