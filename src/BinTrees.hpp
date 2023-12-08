#ifndef BINTREES_HPP
#define BINTREES_HPP

template <typename T>
class BinTreeNode {
public:
    T value;
    BinTreeNode<T> *left, *right, *middle;

    BinTreeNode(T value) : value(value), left(nullptr), right(nullptr), middle(nullptr) {}

    // sum values of the node and its children
    T sumValue() const {
        T sum = value;
        if (left) sum += left->value;
        if (middle) sum += middle->value;
        if (right) sum += right->value;
        return sum;
    }
};

template <typename T>
class BinTree {
public:

    BinTree() : root(nullptr) {}
    ~BinTree() {
        deleteSubTree(root);
    }

    // sum values of the tree
    T sumValue() const {
        if (root) return root->sumValue();
        else return T();
    }

private:
    BinTreeNode<T> *root;
    
    // set sums for layer 0 of the tree
    void set_bintree_sums_layer0(std::vector<T>& sums) {
        if (!root) return;

        sums.clear();
        sums.push_back(root->sumValue()); // Sum values of the root and its children

        // If you have further layers, you would continue the traversal here.
    }

    void deleteSubTree(BinTreeNode<T>* node) {
        if (node) {
            deleteSubTree(node->left);
            deleteSubTree(node->middle);
            deleteSubTree(node->right);
            delete node;
        }
    }

};

#endif // BINTREES_HPP