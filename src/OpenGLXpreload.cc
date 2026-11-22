#include  <GL/glx.h>
#include  <GL/gl.h>
#include  <unistd.h>

static int attributeList[] = { GLX_RGBA, None };

void OpenGLXpreload() {
    Display *dpy;
    XVisualInfo *vi;
    /* get a connection */
    dpy = XOpenDisplay(0);
    /* get an appropriate visual */
    vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeList);
    (void)vi;
}
