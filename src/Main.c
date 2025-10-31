#if defined(__linux__) && !defined(_WIN32)
    #include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#elif defined(_WIN32) || defined(_WIN64)
    #include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
    //#include "F:/home/codeleaded/System/Static/Library/OMML.h"
    //#include "F:/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#elif defined(__APPLE__)
    #error "Apple not supported!"
#else
    #error "Platform not supported!"
#endif



void Setup(AlxWindow* w){
    
}
void Update(AlxWindow* w){
    

    Clear(BLACK);

    
}
void Delete(AlxWindow* w){
    
}
int main(){
    if(Create("Font Renderer",1900,1000,1,1,Setup,Update,Delete))
       Start();
    return 0;
}