#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>
#include <memory>

struct Node {
    std::shared_ptr<uint8_t> value_;
    std::shared_ptr<Node> left_, right_;
};

// HuffmanTree decoder for DHT section.
class HuffmanTree {

public:
    // code_lengths is the array of size no more than 16 with number of
    // terminated nodes in the Huffman tree.
    // values are the values of the terminated nodes in the consecutive
    // level order.

    void Build(const std::vector<uint8_t>& code_lengths, const std::vector<uint8_t>& values);

    // Moves the state of the huffman tree by |bit|. If the node is terminated,
    // returns true and overwrites |value|. If it is intermediate, returns false
    // and value is unmodified.
    bool Move(bool bit, int& value);

    void Clean(std::shared_ptr<Node>&);

    size_t size = 0;

private:
    std::shared_ptr<Node> root_;
    std::shared_ptr<Node> curr_state_;
};
