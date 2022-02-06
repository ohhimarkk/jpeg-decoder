#include "huffman.h"
#include <exception>
#include <iostream>

void HuffmanTree::Clean(std::shared_ptr<Node> &curr_node) {
    if (curr_node == nullptr) {
        return;
    }
    Clean(curr_node->left_);
    Clean(curr_node->right_);
    curr_node.reset();
}

void HuffmanTree::Build(const std::vector<uint8_t> &code_lengths,
                        const std::vector<uint8_t> &values) {
    if (code_lengths.size() > 16) {
        throw std::invalid_argument("Invalid levels amount");
    }
    std::vector<std::shared_ptr<Node>> curr_level, next_level;
    Clean(root_);
    curr_level.push_back(root_ = std::make_shared<Node>());
    curr_state_ = root_;
    size_t free = 0, curr_value = 0;
    size_t last_level = code_lengths.size();
    if (!code_lengths.empty()) {
        while (last_level && code_lengths[--last_level] == 0) {
        }
        ++last_level;
    }
    for (size_t level = 0; level < last_level; ++level) {
        next_level.reserve(curr_level.size() * 2);
        for (size_t i = free; i != curr_level.size(); ++i) {
            next_level.push_back(curr_level[i]->left_ = std::make_shared<Node>());
            next_level.push_back(curr_level[i]->right_ = std::make_shared<Node>());
        }
        curr_level = std::move(next_level);
        next_level.clear();

        if (curr_level.size() < code_lengths[level]) {
            throw std::invalid_argument("Invalid code length");
        }

        for (uint8_t i = 0; i != code_lengths[level]; ++i) {
            if (curr_value >= values.size()) {
                throw std::invalid_argument("Invalid values size");
            }
            curr_level[i]->value_ = std::make_shared<uint8_t>(values[curr_value++]);
        }
        free = code_lengths[level];
    }
}

bool HuffmanTree::Move(bool bit, int &value) {
    if (curr_state_ == nullptr) {
        throw std::invalid_argument("Invalid move");
    }
    if (bit) {
        curr_state_ = curr_state_->right_;
    } else {
        curr_state_ = curr_state_->left_;
    }
    if (curr_state_ == nullptr) {
        throw std::invalid_argument("Invalid move");
    } else if (curr_state_->value_ != nullptr) {
        value = *(curr_state_->value_);
        curr_state_ = root_;
        return true;
    }
    return false;
}
