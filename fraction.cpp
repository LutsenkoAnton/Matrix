#include "fraction.h"
Fraction::Fraction() : numerator_(0), denominator_(1) {}
Fraction::Fraction(int64_t a) : numerator_(a), denominator_(1) {}
Fraction::Fraction(int64_t a, int64_t b) : numerator_(a), denominator_(b) {
    Shorten();
}

std::strong_ordering Fraction::operator<=>(const Fraction& other) const {
    return numerator_ * other.denominator_ <=> denominator_ * other.numerator_;
}

Fraction Fraction::operator+(const Fraction& other) const {
    return Fraction(numerator_ * other.denominator_ + other.numerator_ * denominator_,
                    denominator_ * other.denominator_);
}

Fraction Fraction::operator-(const Fraction& other) const {
    return Fraction(numerator_ * other.denominator_ - other.numerator_ * denominator_,
                    denominator_ * other.denominator_);
}

Fraction Fraction::operator-() const {
    return Fraction(-numerator_, denominator_);
}

Fraction Fraction::operator*(const Fraction& other) const {
    return Fraction(numerator_ * other.numerator_, denominator_ * other.denominator_);
}

Fraction Fraction::operator/(const Fraction& other) const {
    return Fraction(numerator_ * other.denominator_, denominator_ * other.numerator_);
}

Fraction& Fraction::operator+=(const Fraction& other) {
    (*this) = (*this) + other;
    return *this;
}

Fraction& Fraction::operator-=(const Fraction& other) {
    return *this += -other;
}

Fraction& Fraction::operator*=(const Fraction& other) {
    numerator_ *= other.numerator_;
    denominator_ *= other.denominator_;
    Shorten();
    return *this;
}

Fraction& Fraction::operator/=(const Fraction& other) {
    numerator_ *= other.denominator_;
    denominator_ *= other.numerator_;
    Shorten();
    return *this;
}

void Fraction::Shorten() {
    int64_t d = std::gcd(numerator_, denominator_);
    numerator_ /= d;
    denominator_ /= d;
    if (denominator_ < 0) {
        numerator_ *= -1;
        denominator_ *= -1;
    }
}

std::ostream& operator<<(std::ostream& stream, const Fraction& fraction) {
    if (fraction.denominator_ == 1) {
        stream << fraction.numerator_;
        return stream;
    }
    stream << fraction.numerator_ << '/' << fraction.denominator_;
    return stream;
}
