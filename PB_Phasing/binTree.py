#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

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
            raise TypeError('p must be proper Position type.')
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
        '''Creat a new left child for Position p, storing element e.
        Return the Position of new node.
        Raise ValueError if Position p is invalid or p already has a left child.
        '''
        node = self.__validate(p)
        if node.getLeft() is not None:
            raise ValueError('Left child exists.')
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
            for other in self.__sutree_preorder(c):
                yield other

    def __subtree_postorder(self, p):
        '''Generate a postorder iteration of positions in subtree rooted at P.'''
        for c in self.children(p):
            for other in self.__subtree_postorder(c):
                yield other
        yield p

    def __subtree_inorder(self, p):
        '''Generate a preorder iteration of positions in subtree rooted at p.(only for binary tree.'''
        if self.left(p) is not None:
            for other in self.__subtree_inorder(self.left(p)):
                yield other
        yield p
        if self.right(p) is not None:
            for other in self.__subtree_inorder(self.right(p)):
                yield other

    def inorder_indent(self, p, d):
        '''Print inorder representation of subtree of T rooted at p at depth d.'''
        if not self.is_empty():
            for other in self.__subtree_inorder(p):
                dep = self.depth(other) + d
                print(2*dep*' ', '--', other.getElement())

    def delete_subtree(self, p):
        '''Delete the node at Position p, and its child, if any.
        Only work in post order.
        '''
        if not self.is_empty():
            for other in self.__subtree_postorder(p):
                node = other.getNode()
                node.setElement(None)
                node.setLeft(None)
                node.setRight(None)
                node.setParent(None)

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

    def depth(self, p):
        '''Return the number of levels separation Position p from the root.'''
        if self.is_root(p):
            return 0
        else:
            return 1 + self.depth(self.parent(p))

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

    def __height(self, p=None):
        '''Return the height of the subtree rooted at Position p.
        If p is None return the height of entire tree.
        '''
        if p is None:
            p = self.root()
        if self.is_leaf(p):
            return 0
        return 1 + max(self. height(c) for c in self.children(p))

if __name__ == '__main__':
    import cpuCheck, os
    pid = os.getpid()
    print('origin:', cpuCheck.mem(pid))
    '''
    Test for binary tree init and print out
    '''
    lbt = LinkedBinaryTree()
    p0 = lbt.add_root('/')
    p00 = lbt.add_left(p0,'Document')
    p000 = lbt.add_left(p00,'PDFs')
    p001 = lbt.add_right(p00,'PNGs')
    p01 = lbt.add_right(p0,'Desktop')
    p010 = lbt.add_left(p01,'Wallpapers')
    p011 = lbt.add_right(p01,'Folds')
    for x in lbt:
        print(x)
    lbt.inorder_indent(lbt.root(), 0)
    '''
    Test for memory realse
    '''
    lbt = LinkedBinaryTree()
    p = lbt.add_root('/')
    for x in range(0,100):
        p = lbt.add_left(p, [0]*1000)
    p = lbt.root()
    for x in range(0,100):
        p = lbt.add_right(p, [1]*1000)
    print('init tree:', cpuCheck.mem(pid))
    lbt.delete_subtree(lbt.root())
    print('delete tree:', cpuCheck.mem(pid))
