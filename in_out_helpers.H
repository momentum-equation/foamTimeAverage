template<class FieldType>
bool headerChecker
(
    IOobject& fieldHeader
)
{
    bool check = false;
    check = fieldHeader.typeHeaderOk<FieldType>(true);
    return check;
}

IOobject generateDictionary
( 
    word& dictName,
    fvMesh& mesh   
)
{
    IOobject dictIO
    (
        dictName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ
    );

    return dictIO;
}  
