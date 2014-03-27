#ifndef TRIPLET_H
#define TRIPLET_H

template <class T>
class Triplet
{
public:
    Triplet(T first, T second, T third):first(first),second(second),third(third){}
    Triplet(){}

 T first;
 T second;
 T third;
};

#endif // TRIPLET_H
