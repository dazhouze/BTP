#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Position class: 
    __init__(): 
        __container: a ref to a instance of the Positional list   _      _      _      _ 
                                                               <-|_|-><-|_|-><-|_|-><-|_|->   
        __node: a ref to a intance of __Node class    _    
                                                   <-|_|->  
    __eq__() __ne__()                      
    element()
    __make_position()
    __validate()

__Node class: is the base unit.   
       _
    <-|_|->
         __init__()
         getPrev()    setPrev()
         getNext()    setNext()
         getElement() setElement()

PositionalList class: positional deque.
       _      _      _      _ 
    <-|_|-><-|_|-><-|_|-><-|_|->
    basic func:
        __init__()
        __len__()
        __iter__()
        is_empty()
        first()    last()
        before(p)  after(p)
    insert + delete func:
        __insert_between()   __delete_between()
        add_first()          add_last()
        add_before()         add_after()
        delet_()
        replace()
'''
class PositionalList(object):
    '''A sequential container of elements allowing positional access.'''

    ##### Position class#####
    class Position(object):
        '''An abstraction representing the location of a single element.'''
        def __init__(self, container, node):
            '''Constructor should not be invoked by user.'''
            self.__container = container # instance of PositionList class
            self.__node = node # instance of __Node class
            
        def getContainer(self):
            return self.__container

        def getNode(self):
            return self.__node

        def element(self):
            '''Return the element stored at this Position.'''
            return self.getNode().getElement()

        def __eq__(self, other):
            '''Return True if other is a Position represeting the same location.'''
            return type(other) is type(self) and other.getNode() is self.getNode()

        def __ne__(self, other):
            '''Retrun True if other does not represent the same loaction.'''
            return not (self == other)

    ##### utility method #####
    def __validate(self, p):
        '''Return position's node, or raise approprate error if invalid.'''
        if not isinstance(p, self.Position):
            raise TypeError('p must be proper Position type')
        if p.getContainer() is not self:
            raise ValueError('p does not belong to this container')
        if p.getNode().getNext() is None:
            raise ValueError('p is no longer valid')
        return p.getNode()

    def __make_position(self, node):
        '''Return Position instance for given node (or None if sentinel).'''
        if node is self.__header or node is self.__trailer:
            return None
        return self.Position(self, node)

    ##### __Node class #####
    class __Node(object):
        '''Lightweigth, nonpublic class for storing a double linked node.'''
        __slots__ = '__element', '__prev', '__next'

        def __init__(self, e, p, n):
            self.__element = e
            self.__prev = p
            self.__next = n

        def getPrev(self):
            return self.__prev

        def getNext(self):
            return self.__next

        def getElement(self):
            return self.__element

        def setPrev(self, p):
            self.__prev = p

        def setNext(self, n):
            self.__next = n

        def setElement(self, e):
            self.__element = e

    ##### Positional list class #####
    def __init__(self):
        '''Creat an empty list'''
        self.__header = self.__Node(None, None, None)
        self.__trailer = self.__Node(None, None, None)
        self.__header.setNext(self.__trailer)
        self.__trailer.setPrev(self.__header)
        self.__size = 0

    def __len__(self):
        '''Return the number of elements in the list.'''
        return self.__size

    def is_empty(self):
        '''Return True if the list is empty.'''
        return self.__size == 0


    ##### accessors #####
    def first(self):
        '''Return the first Position in the list (or None if list is empty).'''
        return self.__make_position(self.__header.getNext())

    def last(self):
        '''Return the first Position in the list (or None if list is empty).'''
        return self.__make_position(self.__trailer.getPrev())

    def before(self, p):
        '''Return the Position just before Position p (or None if p is first).'''
        node = self.__validate(p)
        return self.__make_position(node.getPrev())
    
    def after(self, p):
        '''Return the Position just after Position p (or None if p is last).'''
        node = self.__validate(p)
        return self.__make_position(node.getNext())

    def __iter__(self):
        '''Generatea forward iteration of the elements of the list.'''
        cursor = self.first()
        while cursor is not None:
            yield cursor.element()
            cursor = self.after(cursor)

    ##### mutators #####
    def __insert_between(self, e, predecessor, successor):
        '''Add element e between two existing nodes and return new node.'''
        newest = self.__Node(e, predecessor, successor)
        predecessor.setNext(newest)
        successor.setPrev(newest)
        self.__size += 1
        return self.__make_position(newest)

    def __delete_node(self, node):
        '''Delete nonsentinel node from the list and returen its element.'''
        predecessor = node.getPrev()
        successor = node.getNext()
        predecessor.setNext(successor)
        successor.setPrev(predecessor)
        self.__size -= 1
        element = node.getElement()
        node.setPrev(None)
        node.setNext(None)
        node.setElement(None)
        return element

    def add_first(self, e):
        '''Insert element e at the font  of the list and return new Postion.'''
        return self.__insert_between(e, self.__header, self.__header.getNext())

    def add_last(self, e):
        '''Insert element e at the back of the list and return new position.'''
        return self.__insert_between(e, self.__trailer.getPrev(), self.__trailer)
    
    def add_before(self, p, e):
        '''Insert element e into list after Positon p and return new Postion.'''
        original = self.__validate(p)
        return self.__insert_between(e, original.getPrev(), original)

    def add_after(self, p, e):
        '''Insert element e into list after Position pand return new Position.'''
        original = self.__validate(p)
        return self.__insert_between(e, original, original.getNext())

    def delete(self, p):
        '''Remove and return the elemet at Position p.'''
        original = self.__validate(p)
        return self.__delete_node(original)

    def replace(self, p, e):
        '''
        Replase the element at Position p.
        Retrun the element formerly at Position p.
        '''
        original = self.__validate(p)
        old_value = orginal.getElement()
        original.setElement(e)
        return old_value

if __name__ == '__main__':
    PL = PositionalList()
    p1 = PL.add_first('H')
    p2 = PL.add_after(p1, 'E')
    p5 = PL.add_last('O')
    p4 = PL.add_before(p5, 'L')
    p3 = PL.add_before(p4, 'L')
    p = PL.last()
    for x in ' WORLD':
        p = PL.add_after(p, x)
    print('length of PositionalList:%d'%len(PL))
    print(PL.first().getNode().getElement())
    for x in PL:
        print(x, end = ' ')
    print('')
