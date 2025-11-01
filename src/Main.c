#if defined(__linux__) && !defined(_WIN32)
    #include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
    #include "/home/codeleaded/System/Static/Library/StdFont.h"
#elif defined(_WIN32) || defined(_WIN64)
    #include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
    #include "/home/codeleaded/System/Static/Library/StdFont.h"
    //#include "F:/home/codeleaded/System/Static/Library/OMML.h"
    //#include "F:/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#elif defined(__APPLE__)
    #error "Apple not supported!"
#else
    #error "Platform not supported!"
#endif



StdFont stdfont;

void Setup(AlxWindow* w){
    stdfont = StdFont_New();
}
void Update(AlxWindow* w){
    

    Clear(BLACK);

    StdFont_Render_CStr(WINDOW_STD_ARGS,&stdfont,"Hello World",0.0f,0.0f,WHITE);
}
void Delete(AlxWindow* w){
    StdFont_Free(&stdfont);
}
int main(){
    if(Create("Font Renderer",1900,1000,1,1,Setup,Update,Delete))
       Start();
    return 0;
}