#ifndef G4REFERENCECAST_HH
#define G4REFERENCECAST_HH

template<typename ReturnType, typename OriginalType>
ReturnType&
reference_cast(OriginalType& source)
{
    return (ReturnType &) source;
}

template<typename ReturnType, typename OriginalType>
ReturnType&
reference_cast(ReturnType&, OriginalType& source)
{
    return (ReturnType &) source;
}

#endif // G4REFERENCECAST_HH
