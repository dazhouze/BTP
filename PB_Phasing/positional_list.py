#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Double linked list for reads information.'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Position class:
    __init__():
        __container: a ref to a instance of the Positional list
        __node: a ref to a intance of _Node class    p
                                                   <-|_|->
    __eq__() __ne__()
    element()
    __make_position()
    __validate()

_Node class: is the base unit.
       _
    <-|_|->
         __init__()
         get_prev()    set_prev()
         get_next()    set_next()
         get_element() set_element()

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
            self.__node = node # instance of _Node class

        def get_container(self):
            ''' Return container'''
            return self.__container

        def get_node(self):
            ''' Return node'''
            return self.__node

        def get_element(self):
            '''Return the element stored at this Position.'''
            return self.get_node().get_element()

        def __eq__(self, other):
            '''Return True if other is a Position represeting the same location.'''
            return type(other) is type(self) and other.get_node() is self.get_node()

        def __ne__(self, other):
            '''Retrun True if other does not represent the same loaction.'''
            return not self == other

    ##### utility method #####
    def __validate(self, p):
        '''Return position's node, or raise approprate error if invalid.'''
        if not isinstance(p, self.Position):
            raise TypeError('p must be proper Position type')
        if p.get_container() is not self:
            raise ValueError('p does not belong to this container')
        if p.get_node().get_next() is None:
            raise ValueError('p is no longer valid')
        return p.get_node()

    def __make_position(self, node):
        '''Return Position instance for given node (or None if sentinel).'''
        if node is self.__header or node is self.__trailer:
            return None
        return self.Position(self, node)

    ##### _Node class #####
    class _Node(object):
        '''Lightweigth, nonpublic class for storing a double linked node.'''
        __slots__ = '__element', '__prev', '__next'

        def __init__(self, e, p, n):
            self.__element = e
            self.__prev = p
            self.__next = n

        def get_prev(self):
            ''' Return previous node Position.'''
            return self.__prev

        def get_next(self):
            ''' Return next node Position.'''
            return self.__next

        def get_element(self):
            ''' Return element Position.'''
            return self.__element

        def set_prev(self, p):
            ''' Set previous node Position.'''
            self.__prev = p

        def set_next(self, n):
            ''' Set next node Position.'''
            self.__next = n

        def set_element(self, e):
            ''' Set element Position.'''
            self.__element = e

    ##### Positional list class #####
    def __init__(self):
        '''Creat an empty list'''
        self.__header = self._Node(None, None, None)
        self.__trailer = self._Node(None, None, None)
        self.__header.set_next(self.__trailer)
        self.__trailer.set_prev(self.__header)
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
        return self.__make_position(self.__header.get_next())

    def last(self):
        '''Return the first Position in the list (or None if list is empty).'''
        return self.__make_position(self.__trailer.get_prev())

    def before(self, p):
        '''Return the Position just before Position p (or None if p is first).'''
        node = self.__validate(p)
        return self.__make_position(node.get_prev())

    def after(self, p):
        '''Return the Position just after Position p (or None if p is last).'''
        node = self.__validate(p)
        return self.__make_position(node.get_next())

    def __iter__(self):
        '''Generatea forward iteration of the elements of the list.'''
        cursor = self.first()
        while cursor is not None:
            yield cursor.get_element()
            cursor = self.after(cursor)

    ##### mutators #####
    def __insert_between(self, e, predecessor, successor):
        '''Add element e between two existing nodes and return new node.'''
        newest = self._Node(e, predecessor, successor)
        predecessor.set_next(newest)
        successor.set_prev(newest)
        self.__size += 1
        return self.__make_position(newest)

    def __delete_node(self, node):
        '''Delete nonsentinel node from the list and returen its element.'''
        predecessor = node.get_prev()
        successor = node.get_next()
        predecessor.set_next(successor)
        successor.set_prev(predecessor)
        self.__size -= 1
        element = node.get_element()
        node.set_prev(None)
        node.set_next(None)
        node.set_element(None)
        return element

    def add_first(self, e):
        '''Insert element e at the font  of the list and return new Postion.'''
        return self.__insert_between(e, self.__header, self.__header.get_next())

    def add_last(self, e):
        '''Insert element e at the back of the list and return new position.'''
        return self.__insert_between(e, self.__trailer.get_prev(), self.__trailer)

    def add_before(self, p, e):
        '''Insert element e into list after Positon p and return new Postion.'''
        original = self.__validate(p)
        return self.__insert_between(e, original.get_prev(), original)

    def add_after(self, p, e):
        '''Insert element e into list after Position pand return new Position.'''
        original = self.__validate(p)
        return self.__insert_between(e, original, original.get_next())

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
        old_value = orginal.get_element()
        original.set_element(e)
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
    print(PL.first().get_node().get_element())
    for x in PL:
        print(x, end=' ')
    print('')
    p = PL.first()
    while p:
        print(p.get_element(), end='')
        p = PL.after(p)
    print('')
    n = PL.delete(p1)
    p = PL.first()
    while p:
        print(p.get_element(), end='')
        p = PL.after(p)
