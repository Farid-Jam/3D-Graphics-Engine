#include <SDL.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

int FOV = 90;
int distance = 20;
int windowWidth = 1280;
int windowHeight = 700;
int running = 1;
int culling = 0;
int fillMode = 0;
int nObjects = 2;
double movementSpeed = 5;
double camYaw = 0;
double camPitch = 0;
char* bunny = "models/bunny.obj";
char* head = "models/head.obj";

// Z-buffer for depth testing
float *zBuffer = NULL;

// Structure to store 3-dimensional points
typedef struct point3D
{
    double x;
    double y;
    double z;
} point3D;

// Structure to store 2-dimensional points with depth
typedef struct point2D
{
    double x;
    double y;
    double z; // depth value
} point2D;

// Structure to store a triangle face of an object
typedef struct triangle
{
    point3D p[3];
    point3D normal;
    double dot;
} triangle;

typedef struct object
{
    int nTriangles;
    triangle* triangles;
} object;

// Initialize or resize Z-buffer
void initZBuffer(int width, int height)
{
    if (zBuffer) {
        free(zBuffer);
    }
    zBuffer = (float*)malloc(width * height * sizeof(float));
}

// Clear Z-buffer with far plane value
void clearZBuffer(int width, int height)
{
    for (int i = 0; i < width * height; i++) {
        zBuffer[i] = -FLT_MAX; // Initialize to negative infinity
    }
}

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

// Method to project a 3-dimensional point onto a 2-dimensional plane with depth
point2D project(point3D p, point3D cam)
{
    // Translate to camera space
    double x = p.x - cam.x;
    double y = p.y - cam.y;
    double z = p.z - cam.z;

    // Apply yaw rotation around global Y
    double xz = sqrt(x*x + z*z);
    double theta = atan2(z, x) - camYaw;
    x = xz * cos(theta);
    z = xz * sin(theta);

    // Apply pitch rotation around camera local X
    double yz = sqrt(y*y + z*z);
    theta = atan2(z, y) - camPitch;
    y = yz * cos(theta);
    z = yz * sin(theta);

    double fov_rad = (FOV * M_PI) / 180;
    double factor = 1 / tan(fov_rad / 2);
    if (z == 0) z = 0.001;

    // Return the projected point with depth
    point2D projected;
    projected.x = x / z * factor * windowWidth + windowWidth/2;
    projected.y = y / z * factor * windowWidth + windowHeight/2;
    projected.z = z; // Store depth for z-buffering
    return projected;
}

// Interpolate depth values across scanline
void fillTriangleZBuffer(SDL_Renderer *renderer, point2D p1, point2D p2, point2D p3)
{
    // Sort points by y-coordinate ascending order
    if (p2.y < p1.y) { point2D temp = p1; p1 = p2; p2 = temp; }
    if (p3.y < p1.y) { point2D temp = p1; p1 = p3; p3 = temp; }
    if (p3.y < p2.y) { point2D temp = p2; p2 = p3; p3 = temp; }

    // Determine line slopes
    double dx1 = (p2.y - p1.y) > 0 ? (p2.x - p1.x) / (p2.y - p1.y) : 0;
    double dx2 = (p3.y - p1.y) > 0 ? (p3.x - p1.x) / (p3.y - p1.y) : 0;
    double dx3 = (p3.y - p2.y) > 0 ? (p3.x - p2.x) / (p3.y - p2.y) : 0;

    // Depth slopes
    double dz1 = (p2.y - p1.y) > 0 ? (p2.z - p1.z) / (p2.y - p1.y) : 0;
    double dz2 = (p3.y - p1.y) > 0 ? (p3.z - p1.z) / (p3.y - p1.y) : 0;
    double dz3 = (p3.y - p2.y) > 0 ? (p3.z - p2.z) / (p3.y - p2.y) : 0;

    // Determine start and end x values
    double sx = p1.x;
    double ex = p1.x;
    double sz = p1.z;
    double ez = p1.z;

    // Determine start and end y values
    int sy = (int)ceil(p1.y);
    int ey = (int)ceil(p2.y);

    // Fill from p1 to p2
    for (int y = sy; y < ey; y++) {
        if (y < 0 || y >= windowHeight) continue;
        
        double dy = y - p1.y;
        int startX = (int)(sx + dy * dx1);
        int endX = (int)(ex + dy * dx2);
        float startZ = (float)(sz + dy * dz1);
        float endZ = (float)(ez + dy * dz2);
        
        if (startX > endX) { 
            int tempX = startX; startX = endX; endX = tempX; 
            float tempZ = startZ; startZ = endZ; endZ = tempZ;
        }

        // Draw scanline with z-buffer test
        for (int x = startX; x <= endX; x++) {
            if (x < 0 || x >= windowWidth) continue;
            
            int idx = y * windowWidth + x;
            float t = (endX - startX) > 0 ? (float)(x - startX) / (endX - startX) : 0;
            float z = startZ + t * (endZ - startZ);
            
            if (z > zBuffer[idx]) {
                zBuffer[idx] = z;
                SDL_RenderDrawPoint(renderer, x, y);
            }
        }
    }

    // Redetermine start and end values for filling from p2 to p3
    sx = p2.x;
    sz = p2.z;
    sy = (int)ceil(p2.y);
    ey = (int)ceil(p3.y);

    // Fill from p2 to p3
    for (int y = sy; y < ey; y++) {
        if (y < 0 || y >= windowHeight) continue;
        
        double dy1 = y - p2.y;
        double dy2 = y - p1.y;
        int startX = (int)(sx + dy1 * dx3);
        int endX = (int)(ex + dy2 * dx2);
        float startZ = (float)(sz + dy1 * dz3);
        float endZ = (float)(ez + dy2 * dz2);
        
        if (startX > endX) { 
            int tempX = startX; startX = endX; endX = tempX; 
            float tempZ = startZ; startZ = endZ; endZ = tempZ;
        }

        // Draw scanline with z-buffer test
        for (int x = startX; x <= endX; x++) {
            if (x < 0 || x >= windowWidth) continue;
            
            int idx = y * windowWidth + x;
            float t = (endX - startX) > 0 ? (float)(x - startX) / (endX - startX) : 0;
            float z = startZ + t * (endZ - startZ);
            
            if (z > zBuffer[idx]) {
                zBuffer[idx] = z;
                SDL_RenderDrawPoint(renderer, x, y);
            }
        }
    }
}

// Function to fill cube faces with z-buffering
void fillFace(SDL_Renderer *renderer, triangle tri, point3D cam)
{
    // Project vertexes onto 2-D plane
    point2D v1 = project(tri.p[0], cam);
    point2D v2 = project(tri.p[1], cam);
    point2D v3 = project(tri.p[2], cam);

    // Determine color of face based off relativity to camera
    double brightness = tri.dot;
    if (brightness > 1.0) brightness = 1.0;
    if (brightness < 0.0) brightness = 0.0;
    Uint8 r = (Uint8)(143 * tri.dot);
    Uint8 g = (Uint8)(153 * tri.dot);
    Uint8 b = (Uint8)(251* tri.dot);
    SDL_SetRenderDrawColor(renderer, r, g, b, 255);

    // Fill triangle with z-buffer
    fillTriangleZBuffer(renderer, v1, v2, v3);
}

// .obj data loader
object loadOBJ(const char* filename)
{
    FILE* f = fopen(filename, "r");
    if (!f){
        printf("Failed to open %s\n", filename); 
        object empty = {0, NULL};
        return empty;
    }

    char line[512];
    int vertexCapacity = 1000000;
    int triCapacity = 2000000;
    int nVertices = 0;
    int nTriangles = 0;

    point3D* vert = malloc(sizeof(point3D) * vertexCapacity);
    triangle* tri = malloc(sizeof(triangle) * triCapacity);

    while (fgets(line, sizeof(line), f)){
        if (strncmp(line, "v ", 2) == 0){
            if (nVertices >= vertexCapacity){
                vertexCapacity *= 2;
                vert = realloc(vert, sizeof(point3D) * vertexCapacity);
            }
            point3D v;
            sscanf(line, "v %lf %lf %lf", &v.x, &v.y, &v.z);
            vert[nVertices++] = v;
        } else if (strncmp(line, "f ", 2) == 0){
            int idx[10]; // support polygons up to 10 vertices
            int count = 0;
            char* ptr = line + 2;
            while (*ptr && count < 10){
                int vi;
                if (sscanf(ptr, "%d", &vi) != 1) break;
                idx[count++] = vi;
                char* next = strchr(ptr, ' ');
                if (!next) break;
                ptr = next + 1;
            }

            // triangulate polygon (fan method)
            for (int i = 1; i < count - 1; i++){
                if (nTriangles >= triCapacity){
                    triCapacity *= 2;
                    tri = realloc(tri, sizeof(triangle) * triCapacity);
                }
                tri[nTriangles] = (triangle){
                    vert[idx[0]-1],
                    vert[idx[i]-1],
                    vert[idx[i+1]-1],
                    {0,0,0},
                    0.0
                };
                nTriangles++;
            }
        }
    }

    // free memory
    fclose(f);
    free(vert); 

    object obj;
    obj.nTriangles = nTriangles;
    obj.triangles = tri;
    return obj;
}

// Function to render objects based on depth level
void renderObject(object obj, point3D cam, SDL_Renderer *renderer, int fillMode, int culling) 
{
    for (int i = 0; i < obj.nTriangles; i++){
        obj.triangles[i].normal = findNormal(obj.triangles[i]);
        obj.triangles[i].dot = findDot(obj.triangles[i], cam);

        if (fillMode){
            // In fill mode, use z-buffer for all triangles
            if (!culling || obj.triangles[i].dot > 0) {
                fillFace(renderer, obj.triangles[i], cam);
            }
        } else {
            // Wireframe mode (no z-buffer needed)
            point2D v1 = project(obj.triangles[i].p[0], cam);
            point2D v2 = project(obj.triangles[i].p[1], cam);
            point2D v3 = project(obj.triangles[i].p[2], cam);

            if ((obj.triangles[i].dot > 0 && culling) || !culling){
                SDL_RenderDrawLine(renderer, v1.x, v1.y, v2.x, v2.y);
                SDL_RenderDrawLine(renderer, v2.x, v2.y, v3.x, v3.y);
                SDL_RenderDrawLine(renderer, v3.x, v3.y, v1.x, v1.y);
            }
        }
    }
}


int main()
{
    // Initialize system
    SDL_Init(SDL_INIT_VIDEO);
    SDL_ShowCursor(SDL_DISABLE);
    SDL_SetRelativeMouseMode(SDL_TRUE);

    // Initialize camera position
    point3D camera = {0, 0, 0};

    // Create window
    SDL_Window *window = SDL_CreateWindow("Graphics3D", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
    windowWidth, windowHeight, SDL_WINDOW_SHOWN);
    SDL_SetWindowResizable(window, SDL_TRUE);

    // Grab renderer 
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    // Initialize Z-buffer
    initZBuffer(windowWidth, windowHeight);

    // Load object
    object objects[nObjects];
    objects[0] = loadOBJ(bunny);
    objects[1] = loadOBJ(head);

    // After loading and scaling the modelz, calculate each models center
    for (int obj = 0; obj < nObjects; obj++ ){
        double scale = 50.0;
        point3D center = {0, 0, 0};

        for (int i = 0; i < objects[obj].nTriangles; i++) {
            for (int j = 0; j < 3; j++) {
                objects[obj].triangles[i].p[j].x *= scale;
                objects[obj].triangles[i].p[j].y *= scale;
                objects[obj].triangles[i].p[j].z *= scale;

                center.x += objects[obj].triangles[i].p[j].x;
                center.y += objects[obj].triangles[i].p[j].y;
                center.z += objects[obj].triangles[i].p[j].z;
            }
        }

        // Average to find center
        int totalPoints = objects[obj].nTriangles * 3;
        center.x /= totalPoints;
        center.y /= totalPoints;
        center.z /= totalPoints;

        // Translate model so its center is at origin, then move to distance
        for (int i = 0; i < objects[obj].nTriangles; i++) {
            for (int j = 0; j < 3; j++) {
                objects[obj].triangles[i].p[j].x -= center.x;
                objects[obj].triangles[i].p[j].y -= center.y;
                objects[obj].triangles[i].p[j].z -= center.z;
                objects[obj].triangles[i].p[j].z -= distance; 
            }
        }
    }

    for (int i = 0; i < objects[0].nTriangles; i++) {
            for (int j = 0; j < 3; j++) {
                objects[0].triangles[i].p[j].x -= 16;
                objects[0].triangles[i].p[j].y -= 8;
            }
        }

    // Timing variables for smooth frame-independent movement
    Uint64 lastTime = SDL_GetPerformanceCounter();
    double deltaTime = 0.0;

    // Render loop
    while(running)
    {
        // Calculate delta time
        Uint64 currentTime = SDL_GetPerformanceCounter();
        deltaTime = (double)(currentTime - lastTime) / SDL_GetPerformanceFrequency();
        lastTime = currentTime;

        // Process ALL events in the queue
        SDL_Event event;
        while(SDL_PollEvent(&event))
        {
            // If window is closed, end render loop
            if (event.type == SDL_QUIT)
            {
                running = 0;
            }

            // Toggle fill, culling, and fullscreen modes
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
                if (event.key.keysym.sym == SDLK_ESCAPE)
                {
                    SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
                }
            }

            // If window is resized, re-initialize windowWidth and windowHeight
            if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_RESIZED)
            {
                SDL_GetWindowSize(window, &windowWidth, &windowHeight);
                initZBuffer(windowWidth, windowHeight);
            }

            // Mouse movement for camera rotation
            if (event.type == SDL_MOUSEMOTION)
            {
                camYaw += -1 * event.motion.xrel * 0.002;
                camPitch += -1 * event.motion.yrel * 0.002;
                if (camPitch > M_PI / 2) camPitch = M_PI / 2;
                if (camPitch < -M_PI / 2) camPitch = -M_PI / 2;
            }
        }

        // Get keyboard state for smooth movement
        const Uint8 *keyboard_state_array = SDL_GetKeyboardState(NULL);

        // Calculate forward/right vectors for movement (time-scaled)
        double speed = movementSpeed * deltaTime;
        double forwardX = sin(-1 * camYaw) * speed;
        double forwardZ = cos(-1 * camYaw) * speed;
        double rightX   = cos(-1 * camYaw) * speed;
        double rightZ   = -sin(-1 * camYaw) * speed;

        // Adjust camera position
        if (keyboard_state_array[SDL_SCANCODE_W]) { camera.x -= forwardX; camera.z -= forwardZ; }
        if (keyboard_state_array[SDL_SCANCODE_S]) { camera.x += forwardX; camera.z += forwardZ; }
        if (keyboard_state_array[SDL_SCANCODE_A]) { camera.x += rightX; camera.z += rightZ; }
        if (keyboard_state_array[SDL_SCANCODE_D]) { camera.x -= rightX; camera.z -= rightZ; }
        if (keyboard_state_array[SDL_SCANCODE_SPACE]) camera.y += speed;
        if (keyboard_state_array[SDL_SCANCODE_LSHIFT]) camera.y -= speed;
        
        // Clear Z-buffer for new frame
        clearZBuffer(windowWidth, windowHeight);
        
        // Prepare for next frame by clearing previous render
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // Render every object
        SDL_SetRenderDrawColor(renderer, 204, 204, 255, 225);
        for (int i = 0; i < nObjects; i++) {
            renderObject(objects[i], camera, renderer, fillMode, culling);
        }

        // Display newest render
        SDL_RenderPresent(renderer);
        
        // Cap frame rate to 60 FPS for consistency
        SDL_Delay(1);
    }

    // Free memory
    free(zBuffer);

    // Close everything before shutting down
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}