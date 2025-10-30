#include <SDL.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


int FOV = 90;
int distance = -20;
int windowWidth = 1280;
int windowHeight = 700;
int running = 1;
SDL_Event event;
int numVerts = 0;
int numTriangles = 0;
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
    p->z = newZ + distance;;
}

// Method to rotate three dimensional points about the Y axis 
void rotateY(point3D *p, double theta)
{
    p->z -=distance;
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

// .obj data loader
int loadOBJ(const char* filename, point3D** vertices, triangle** tris, int* nTriangles)
{
    FILE* f = fopen(filename, "r");
    if (!f){
        printf("Failed to open %s\n", filename); 
        return 0;
    }
    char line[256];
    int nVertices = 0;
    *nTriangles = 0;

    point3D* vert = malloc(sizeof(point3D) * 5000);
    triangle* tri = malloc(sizeof(triangle) * 10000);

    while (fgets(line, sizeof(line), f)){
        if (strncmp(line, "v ", 2)==0){
            point3D v;
            sscanf(line, "v %lf %lf %lf", &v.x, &v.y, &v.z);
            vert[nVertices++] = v;
        }else if(strncmp(line,"f ",2)==0){
            int idx[10], count=sscanf(line+2,"%d %d %d %d %d %d %d %d %d %d",&idx[0],&idx[1],&idx[2],&idx[3],&idx[4],&idx[5],&idx[6],&idx[7],&idx[8],&idx[9]);
            for(int i=1;i<count-1;i++){
                tri[*nTriangles]=(triangle){vert[idx[0]-1],vert[idx[i]-1],vert[idx[i+1]-1]};
                (*nTriangles)++;
            }
        }
    }
    fclose(f);
    *vertices = vert;
    *tris = tri;
    *nTriangles = *nTriangles;
    return 1;
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

    // Load object
    point3D* vertices = NULL;
    triangle* triangles = NULL;
    if(!loadOBJ("models/bunny.obj", &vertices, &triangles, &numTriangles)){
        printf("Failed to load OBJ\n"); 
        return 1;
    }

    // After loading and scaling the model, calculate its center
    double scale = 50.0;
    point3D center = {0, 0, 0};

    for (int i = 0; i < numTriangles; i++) {
        for (int j = 0; j < 3; j++) {
            triangles[i].p[j].x *= scale;
            triangles[i].p[j].y *= scale;
            triangles[i].p[j].z *= scale;
            
            center.x += triangles[i].p[j].x;
            center.y += triangles[i].p[j].y;
            center.z += triangles[i].p[j].z;
        }
    }

    // Average to find center
    int totalPoints = numTriangles * 3;
    center.x /= totalPoints;
    center.y /= totalPoints;
    center.z /= totalPoints;

    // Translate model so its center is at origin, then move to distance
    for (int i = 0; i < numTriangles; i++) {
        for (int j = 0; j < 3; j++) {
            triangles[i].p[j].x -= center.x;
            triangles[i].p[j].y -= center.y;
            triangles[i].p[j].z -= center.z;
            triangles[i].p[j].z += distance;  // Move to viewing distance
        }
    }

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

    free(vertices);
    free(triangles);

    // Close everything before shutting down
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}