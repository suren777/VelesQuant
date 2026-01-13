
#ifndef OBJECT_CAHCE_H
#define OBJECT_CAHCE_H

#include <map>
#include <sstream>
#include <iostream>
#include <velesquant/xlw/CellMatrix.h>

namespace velesquant {

#pragma warning (disable:4996)

using namespace xlw;

inline std::string ulToString(unsigned long n)
{
    std::ostringstream theStream;
	theStream << n;
    return theStream.str();
}

inline std::string stripTrailingDot(const std::string &theSource)
{
    size_t posHash = theSource.find('.');
    if(posHash==std::string::npos) 
		return theSource;
    return theSource.substr(0,posHash);
}

inline void check(const std::string &theName)
{
    if(theName.find('.')==std::string::npos)
        throw("Character '.' not allowed in object name");
}


template<typename Object>
class singleton
{
public:
	static Object& instance()
    {
        static Object theObject;
        return theObject;
    }
    virtual ~singleton() {}
protected:
    singleton(){}
private:
    singleton(const singleton& theOther);
    singleton& operator=(const singleton& theOther);
} ;


template<typename T> 
class cachedObject : public singleton<cachedObject<T> >, 
	                 public std::map<std::string, T*> 
{};




class NumberCache : public singleton<NumberCache>,
	                public std::map<std::string, unsigned long>
{};

class ObjectCache : public singleton<ObjectCache>,
	                public std::map<std::string, CellMatrix>
{};

class DiscountCurveCache : public singleton<DiscountCurveCache>,
	                       public std::map<std::string, MyMatrix>
{};

template<typename T>
T* getObject(const std::string &theName)
{
std::string theKey(stripTrailingDot(theName));
	if(cachedObject<T>::instance().find(theKey)==cachedObject<T>::instance().end()) 
		throw("Object "+theName+" not found in Cache");
return cachedObject<T>::instance()[theKey];
};

template<typename T>
std::string placeObject(const std::string &theObjName, T* theObj)
{
std::string theKey(stripTrailingDot(theObjName));
	if (cachedObject<T>::instance().find(theKey)!=cachedObject<T>::instance().end()) 
		cachedObject<T>::instance().erase(cachedObject<T>::instance().find(theKey));

	cachedObject<T>::instance().insert(std::make_pair(theObjName,theObj));
	std::string theName(theObjName);
	theName += ".";
	return theName + (ulToString(NumberCache::instance()[theObjName]++));
};


}
#endif