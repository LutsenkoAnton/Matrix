#pragma once

#include <initializer_list>
#include <iostream>
#include <vector>

class Permutation {
public:
    Permutation() = delete;
    Permutation(size_t n);
    Permutation(const std::vector<size_t>& data);
    Permutation(const std::initializer_list<size_t>& data);
    Permutation(const std::initializer_list<std::initializer_list<size_t>>& cycles);
    Permutation(const Permutation& other) = default;
    Permutation(Permutation&& other) = default;

    Permutation& operator=(const Permutation& other) = default;
    Permutation& operator=(Permutation&& other) = default;
    Permutation& operator=(const std::initializer_list<size_t>& data);
    Permutation& operator=(const std::initializer_list<std::initializer_list<size_t>>& data);

    Permutation operator*(const Permutation& other) const;
    Permutation& operator*=(const Permutation& other);

    bool operator==(const Permutation& other) const = default;
    bool operator!=(const Permutation& other) const = default;

    Permutation Inverse() const;
    Permutation Power(size_t indicator) const;

    size_t operator[](size_t ind) const;

    size_t size() const;

    friend class AllPermutations;
private:
    std::vector<size_t> data;
};

std::ostream& operator<<(std::ostream& stream, const Permutation& p);

class AllPermutations {
public:
    AllPermutations() = delete;
    AllPermutations(size_t n);

    class Iterator {
    public:
        Iterator() = delete;
        Iterator(size_t size, bool is_end);

        Iterator operator++();
        Iterator operator++(int);

        const Permutation& operator*() const;

        bool operator==(const Iterator& other) const;
        bool operator!=(const Iterator& other) const;

    private:
        Permutation p;
        bool is_end = false;
    };

    Iterator begin() const;
    Iterator end() const;

private:
    size_t size;
};
