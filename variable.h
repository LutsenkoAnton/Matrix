#pragma once

#include <compare>
#include <cstddef>

// Variable(0) reserved for constant
struct Variable {
    Variable() {
        id = ++count;
    }
    explicit Variable(size_t id): id(id) {}
    std::strong_ordering operator<=>(const Variable& other) const = default;
    size_t id;
    static size_t count;
};
