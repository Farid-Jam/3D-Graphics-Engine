#include <SDL.h>
#include <math.h>

int FOV = 90;
int distance = 11;
int windowWidth = 1280;
int windowHeight = 700;
int running = 1;
SDL_Event event;
int numTriangles = 12;
int culling = 0; // 0 to turn off, any other number to turn on
int fillMode = 0; // 0 to turn off, any other number to turn on

// Structure to store 3-dimensional points
typedef struct point3D
{
    double x;
    double y;
    double z;
} point3D;

// Structure to store 2-dimensional points
typedef struct point2D
{
    double x;
    double y;
} point2D;

// Structure to store a triangle face of an object
typedef struct triangle
{
    point3D p[3];
    point3D normal;
    double dot;
} triangle;

// Method to normalize vectors
point3D normalize(point3D p)
{
    double norm = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    point3D normalized = {p.x / norm, p.y / norm, p.z / norm};
    return normalized;
}

// Method to determine whether the normal of a face points towards the camera 
point3D findNormal(triangle tri)
{
    point3D edge1 = {tri.p[1].x - tri.p[0].x, tri.p[1].y - tri.p[0].y, tri.p[1].z - tri.p[0].z};
    point3D edge2 = {tri.p[2].x - tri.p[0].x, tri.p[2].y - tri.p[0].y, tri.p[2].z - tri.p[0].z};

    point3D normal = {
        (edge1.y * edge2.z) - (edge1.z * edge2.y),
        (edge1.z * edge2.x) - (edge1.x * edge2.z),
        (edge1.x * edge2.y) - (edge1.y * edge2.x)
    };

    return normalize(normal);
}

// Method to determine how much a faces normal vector differs from a vector from the camera to the center of a face
double findDot(triangle tri, point3D camPos)
{
    point3D center = {
        ((tri.p[0].x + tri.p[1].x + tri.p[2].x) / 3),
        ((tri.p[0].y + tri.p[1].y + tri.p[2].y) / 3), 
        ((tri.p[0].z + tri.p[1].z + tri.p[2].z) / 3)
    };

    point3D vector = {
        camPos.x - center.x,
        camPos.y - center.y,
        camPos.z - center.z
    };

    vector = normalize(vector);

    double dot = vector.x * tri.normal.x + vector.y * tri.normal.y + vector.z * tri.normal.z;

    return dot;
}

// Method to projecta 3-demensional point onto a 2-dimensional plane, returning the 2-dimensional point
point2D project(point3D p)
{
    // Obtain FOV in radians and determine projection factor
    double fov_rad = (FOV * M_PI) / 180;
    double factor = 1 / tan(fov_rad / 2);
    
    // Remove division by zero
    if (p.z == 0) p.z = 0.001;

    // Project points using projection factor and projection formula
    point2D projected;
    projected.x = p.x / p.z * factor * windowWidth + windowWidth / 2;
    projected.y = p.y / p.z * factor * windowWidth + windowHeight / 2;

    // Return the projected point 
    return projected;
}

// Method to rotate three dimensional points about the X axis 
void rotateX(point3D *p, double theta)
{
    p->z -= distance;

    double newY = p->y * cos(theta) - p->z * sin(theta);
    double newZ = p->y * sin(theta) + p->z * cos(theta);

    p->y = newY;
    p->z = newZ + distance;
}

// Method to rotate three dimensional points about the Y axis 
void rotateY(point3D *p, double theta)
{
    p->z -= distance;

    double newX = p->x * cos(theta) + p->z * sin(theta);
    double newZ = -1 * (p->x * sin(theta)) + p->z * cos(theta);

    p->x = newX;
    p->z = newZ + distance;
}

// Function to fill in triangles
void fillTriangle(SDL_Renderer *renderer, point2D p1, point2D p2, point2D p3)
{
    // Sort points by y-coordinate ascending order
    if (p2.y < p1.y) { point2D temp = p1; p1 = p2; p2 = temp; }
    if (p3.y < p1.y) { point2D temp = p1; p1 = p3; p3 = temp; }
    if (p3.y < p2.y) { point2D temp = p2; p2 = p3; p3 = temp; }

    // Determine line slopes
    double dx1 = (p2.y - p1.y) > 0 ? (p2.x - p1.x) / (p2.y - p1.y) : 0;
    double dx2 = (p3.y - p1.y) > 0 ? (p3.x - p1.x) / (p3.y - p1.y) : 0;
    double dx3 = (p3.y - p2.y) > 0 ? (p3.x - p2.x) / (p3.y - p2.y) : 0;

    // Determine start and end x values
    double sx = p1.x;
    double ex = p1.x;

    // Determine start and end y values
    int sy = (int)ceil(p1.y);
    int ey = (int)ceil(p2.y);

    // Fill from p1 to p2
    for (int y = sy; y < ey; y++) {
        int startX = (int)(sx + (y - p1.y) * dx1);
        int endX = (int)(ex + (y - p1.y) * dx2);
        if (startX > endX) { int temp = startX; startX = endX; endX = temp; }
        SDL_RenderDrawLine(renderer, startX, y, endX, y);
    }

    // Redetermine start and end values for filling from p2 to p3
    sx = p2.x;
    sy = (int)ceil(p2.y);
    ey = (int)ceil(p3.y);

    // Fill from p2 to p3
    for (int y = sy; y < ey; y++) {
        int startX = (int)(sx + (y - p2.y) * dx3);
        int endX = (int)(ex + (y - p1.y) * dx2);
        if (startX > endX) { int temp = startX; startX = endX; endX = temp; }
        SDL_RenderDrawLine(renderer, startX, y, endX, y);
    }
}

// Function to fill cube faces
void fillFace(SDL_Renderer *renderer, triangle tri)
{
    // Project vertexes onto 2-D plane
    point2D v1 = project(tri.p[0]);
    point2D v2 = project(tri.p[1]);
    point2D v3 = project(tri.p[2]);

    // Determine color of face based off relativity to camera
    double brightness = tri.dot;
    if (brightness > 1.0) brightness = 1.0;
    if (brightness < 0.0) brightness = 0.0;
    Uint8 r = (Uint8)(143 * tri.dot);
    Uint8 g = (Uint8)(153 * tri.dot);
    Uint8 b = (Uint8)(251* tri.dot);
    SDL_SetRenderDrawColor(renderer, r, g, b, 255);

    // Split face into two triangles and fill them
    fillTriangle(renderer, v1, v2, v3);
}

int main()
{
    // Initialize system
    SDL_Init(SDL_INIT_VIDEO);

    // Initialize camera position
    point3D camera = {0, 0, 0};

    // Create window
    SDL_Window *window = SDL_CreateWindow("Graphics3D", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
    windowWidth, windowHeight, SDL_WINDOW_SHOWN);
    SDL_SetWindowResizable(window, SDL_TRUE);

    // Grab renderer 
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    // Initialize points
    point3D p1 = {-1, -1, -1 + distance};
    point3D p2 = {1, -1, -1 + distance};
    point3D p3 = {1, 1, -1 + distance};
    point3D p4 = {-1, 1, -1 + distance};
    point3D p5 = {-1, -1, 1 + distance};
    point3D p6 = {1, -1, 1 + distance};
    point3D p7 = {1, 1, 1 + distance};
    point3D p8 = {-1, 1, 1 + distance};

    // Initialize array of face struc holding all necessary cube faces
    triangle triangles[12] = {
            {p1, p3, p2},  // CW from outside
        {p1, p4, p3},

        // Front face (+Z)
        {p5, p6, p7},  // CW from outside
        {p5, p7, p8},

        // Bottom face (-Y)
        {p1, p2, p6},  // CW from outside
        {p1, p6, p5},

        // Top face (+Y)
        {p4, p8, p7},  // CW from outside
        {p4, p7, p3},

        // Left face (-X)
        {p1, p8, p4},  // CW from outside
        {p1, p5, p8},

        // Right face (+X)
        {p2, p3, p7},  // CW from outside
        {p2, p7, p6}
    };

    // Render loop
    while(running)
    {
        SDL_PollEvent(&event);

        // If window is closed, end render loop
        if (event.type == SDL_QUIT)
        {
            running = 0;
        }
        
        // Toggle fill and culling modes with 'f' and 'c' keys respectively
        if (event.type == SDL_KEYDOWN)
        {
            if (event.key.keysym.sym == SDLK_f)
            {
                fillMode = (fillMode == 0) ? 1 : 0;
            }
            if (event.key.keysym.sym == SDLK_c)
            {
                culling = (culling == 0) ? 1 : 0;
            }
        }

        if (event.type == SDL_WINDOWEVENT)
        {
            // If window is resized, re-initialize windowWidth and windowHeight
            if (event.window.event == SDL_WINDOWEVENT_RESIZED)
            {
                SDL_GetWindowSize(window, &windowWidth, &windowHeight);
            }
        }

        // If mouse is moved within the window, rotate points in directions x and y accordingly
        if (event.type == SDL_MOUSEMOTION)
        {
            double dx = event.motion.xrel;
            double dy = event.motion.yrel;
            double rotationX = -1 * dx * 2 * M_PI / 1000;
            double rotationY = dy * 2 * M_PI / 1000;
            for (int i = 0; i < numTriangles; i++)
            {
                rotateY(&triangles[i].p[0], rotationX);
                rotateY(&triangles[i].p[1], rotationX);
                rotateY(&triangles[i].p[2], rotationX);

                rotateX(&triangles[i].p[0], rotationY);
                rotateX(&triangles[i].p[1], rotationY);
                rotateX(&triangles[i].p[2], rotationY);

            }
        }
        
        // Prepare for next frame by clearing previous render
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // Draw lines connecting points in every visible face
        SDL_SetRenderDrawColor(renderer, 204, 204, 255, 225);
        for (int i = 0; i < numTriangles; i++)
        {
            point2D v1 = project(triangles[i].p[0]);
            point2D v2 = project(triangles[i].p[1]);
            point2D v3 = project(triangles[i].p[2]);
            triangles[i].normal = findNormal(triangles[i]);
            triangles[i].dot = findDot(triangles[i], camera);
            if (triangles[i].dot > 0 && fillMode){
                fillFace(renderer, triangles[i]);
            } else if (triangles[i].dot > 0 && culling && !fillMode){
                SDL_RenderDrawLine(renderer, v1.x, v1.y, v2.x, v2.y);
                SDL_RenderDrawLine(renderer, v2.x, v2.y, v3.x, v3.y);
                SDL_RenderDrawLine(renderer, v3.x, v3.y, v1.x, v1.y);
            } else if (!culling && !fillMode){
                SDL_RenderDrawLine(renderer, v1.x, v1.y, v2.x, v2.y);
                SDL_RenderDrawLine(renderer, v2.x, v2.y, v3.x, v3.y);
                SDL_RenderDrawLine(renderer, v3.x, v3.y, v1.x, v1.y);
            }
        }

        // Display newest render
        SDL_RenderPresent(renderer);
        SDL_Delay(10);
    }

    // Close everything before shutting down
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}