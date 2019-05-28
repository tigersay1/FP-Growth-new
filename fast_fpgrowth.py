"""
Usage of the module:
# global_set should be list of lists  

from fast_fpgrowth import mine_frequent_itemsets, assocRule
for itemset, support in mine_frequent_itemsets(global_set, min_support, True):
    sets[tuple(itemset)] = support
rules = assocRule(sets, min_confidence)
"""

from collections import defaultdict, namedtuple
import itertools


def mine_frequent_itemsets(transactions, min_support, include_support=False):    
    
    
    def delete_rare_in_transaction(transaction):
        transaction = filter(lambda v: v in items, transaction)
        transaction = sorted(transaction, key=lambda v: items[v], reverse=True)
        return transaction
    
    
    def find_new_set(tree, min_support, suffix):
        for item, nodes in tree.items():
            support = sum(n.count for n in nodes)
            if support >= min_support and item not in suffix:         
                new_set = [item] + suffix
                yield (new_set, support) if include_support else new_set

                # Recursive FP-Tree building
                conditional_tree = conditional_tree_from_paths(tree.prefix_paths(item))
                for s in find_new_set(conditional_tree, min_support, new_set):
                    yield s 
    
    items = defaultdict(lambda: 0) 

    # supporting both types of inputting min_support
    if min_support <= 1:
        min_support = int(min_support * len(transactions))
    
    Tree = FPTree()
    for transaction in transactions:
        for item in transaction:
            items[item] += 1

    # Remove infrequent items
    items = dict((item, support) for item, support in items.items() if support >= min_support)

    # Build FP-tree
    for trans in transactions:
        Tree.add(delete_rare_in_transaction(trans))

    # Search for frequent itemsets
    for itemset in find_new_set(Tree, min_support, []):
        yield itemset

class FPTree(object):
    Route = namedtuple('Route', 'head tail')

    def __init__(self):
        self._root = FPNode(self, None, None)
        self._routes = {}
        

    @property
    def root(self):
        return self._root
    

    def add(self, transaction):
        # Add a transaction to the tree
        point = self._root
        for item in transaction:
            next_point = point.search(item)
            if next_point:
                next_point._count += 1
            else:
                next_point = FPNode(self, item)
                point.make_child(next_point)
                self._make_new_route(next_point)
            point = next_point
            

    def _make_new_route(self, point):
        # Add the given node to the route through all nodes for its item
        try:
            route = self._routes[point.item]
            # add to the tail
            route[1].neighbor = point 
            self._routes[point.item] = self.Route(route[0], point)
        except KeyError:
            # start a new route
            self._routes[point.item] = self.Route(point, point)

            
    def items(self):
        """
        Generate one 2-tuples for each item represented in the tree. The first
        element of the tuple is the item itself, and the second element is a
        generator that will yield the nodes in the tree that belong to the item.
        """
        for item in self._routes.keys():
            yield (item, self.nodes(item))

            
    def nodes(self, item):
        try:
            node = self._routes[item][0]
        except KeyError:
            return

        while node:
            yield node
            node = node.neighbor

            
    def prefix_paths(self, item):
        def collect_path(node):
            path = []
            while node:
                path.append(node)
                node = node.parent
                if node.root:
                    break
            path.reverse()
            return path
        return (collect_path(node) for node in self.nodes(item))


def conditional_tree_from_paths(paths):
    tree = FPTree()
    condition_item = None
    items = set()
    for path in paths:
        if condition_item is None:
            condition_item = path[-1].item
        point = tree.root
        for node in path:
            next_point = point.search(node.item)
            if not next_point:
                # Add a new node to the tree
                items.add(node.item)
                if node.item == condition_item:
                    next_point = FPNode(tree, node.item, node.count)
                else:
                    next_point = FPNode(tree, node.item, 0)
                point.make_child(next_point)
                tree._make_new_route(next_point)
            point = next_point

    # compute the counts 
    for path in tree.prefix_paths(condition_item):
        for node in reversed(path[:-1]):
            node._count += path[-1].count

    return tree

class FPNode(object):

    def __init__(self, tree, item, count=1):
        self._tree = tree
        self._item = item
        self._count = count
        self._parent = None
        self._children = {}
        self._neighbor = None

        
    def make_child(self, child):
        if not child.item in self._children:
            self._children[child.item] = child
            child.parent = self

            
    def search(self, item):
        # check if there is only one child
        try:
            return self._children[item]
        except KeyError:
            return None

    @property
    def tree(self):
        return self._tree

    @property
    def item(self):
        return self._item

    @property
    def count(self):
        return self._count

    @property
    def root(self):
        return self._item is None and self._count is None

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = value

    @property
    def neighbor(self):
        return self._neighbor

    @neighbor.setter
    def neighbor(self, value):
        self._neighbor = value

    @property
    def children(self):
        return tuple(self._children.itervalues())


def subs(l):
    if len(l) == 1:
        return [l]
    x = subs(l[1:])
    return x + [[l[0]] + y for y in x]


def assocRule(freq, min_conf = 0.6):
    assert type(freq) is dict
    result = []
    for item, supp in freq.items():
        for subitem in subs(list(item)):
            if subitem == []:
                continue
            sb = [x for x in item if x not in subitem]
            if sb == []:
                continue
            conf = supp/freq[tuple(subitem)]
            if conf >= min_conf:
                result.append({'from':subitem, 'to':sb, 'supp':supp, 'conf':conf})
    return result