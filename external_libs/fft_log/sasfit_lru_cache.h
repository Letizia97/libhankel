#pragma once
// Fixed-capacity least-recently-used cache, header-only C++17.
//
// Key  — a small POD struct compared byte-for-byte with std::memcmp.
//         Zero-initialise before filling fields (Key k{};) so padding bytes
//         are deterministic and two logically equal keys always compare equal.
// T    — the cached value.  Must be default-constructible; T() is used to
//         reset a slot before the caller fills it in via insert().  Destructor
//         runs automatically when the slot is overwritten or the cache is
//         destroyed — ordinary RAII, no separate release step needed.
//
// Not thread-safe.

#include <cstring>
#include <cstddef>
#include <vector>

template <typename Key, typename T>
class SasfitLruCache {
public:
    explicit SasfitLruCache(int capacity)
        : slots_(static_cast<std::size_t>(capacity > 0 ? capacity : 1)) {}

    // Returns a pointer to the cached value for key, or nullptr on a miss.
    // Marks the entry as most-recently-used on a hit.
    T* lookup(const Key& key) {
        for (auto& s : slots_) {
            if (s.occupied && std::memcmp(&s.key, &key, sizeof(Key)) == 0) {
                s.timestamp = ++clock_;
                return &s.value;
            }
        }
        return nullptr;
    }

    // Returns a reference to a slot ready to be filled for key.
    // Reuses the existing slot if key is already present; otherwise claims an
    // empty slot, or evicts the least-recently-used entry when the cache is
    // full.  The slot's value is reset to T() before this returns.
    T& insert(const Key& key) {
        for (auto& s : slots_) {
            if (s.occupied && std::memcmp(&s.key, &key, sizeof(Key)) == 0) {
                s.value = T();
                s.timestamp = ++clock_;
                return s.value;
            }
        }
        for (auto& s : slots_) {
            if (!s.occupied) {
                s.key       = key;
                s.value     = T();
                s.timestamp = ++clock_;
                s.occupied  = true;
                return s.value;
            }
        }
        Slot* lru = &slots_[0];
        for (auto& s : slots_)
            if (s.timestamp < lru->timestamp) lru = &s;
        lru->key       = key;
        lru->value     = T();
        lru->timestamp = ++clock_;
        return lru->value;
    }

private:
    struct Slot {
        bool          occupied  = false;
        Key           key{};
        T             value{};
        unsigned long timestamp = 0;
    };

    std::vector<Slot> slots_;
    unsigned long     clock_ = 0;
};
