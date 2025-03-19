#include "SDL2/SDL.h"
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <string.h>

// gcc -std=c17 main.c -IC:\Users\frase\OneDrive\Desktop\C\SDL2 -LC:\Users\frase\OneDrive\Desktop\C\SDL2\lib -Wall -lmingw32 -lSDL2main -lSDL2 -o main
// gcc -std=c17 main.c -IC:\Users\frase\OneDrive\Desktop\C\SDL2 -LC:\Users\frase\OneDrive\Desktop\C\SDL2\lib -lgdi32 -lmingw32 -lSDL2main -lSDL2 -o main
// gcc -std=c17 main.c -O3 -flto -march=native -funroll-loops -IC:\Users\frase\OneDrive\Desktop\C\SDL2 -LC:\Users\frase\OneDrive\Desktop\C\SDL2\lib -Wall -lmingw32 -lSDL2main -lSDL2 -o main

typedef struct {
    int x0;
    int x1;
    int y0;
    int y1;
    int is_portal;
    float portal_top;
    float portal_bottom;
    int portal_link;
} wall;
typedef struct {
    float elevation;
    float height;
    size_t count;
    wall walls[];
} sector;
typedef struct {
    size_t count;
    sector *sectors[];
} Level;

typedef struct {
    float x;
    float y;
    float z;
    float h;
    float yaw;
    float focalLength;
    float fov;
} player;

typedef struct {
    float x0;
    float x1;
    float y0;
    float y1;
    float y2;
    float y3;
    float dxCull;
} portalCull;

typedef struct {
    int sector_link;
    portalCull portalBounds;
    int clipped;
} portalRender;

#define PI 3.14159265359
#define DEG2RAD(_d) ((_d) * (PI / 180.0f))
#define RAD2DEG(_d) ((_d) * (180.0f / PI))

void *alloc_dynArray(size_t amnt, size_t target_size, size_t child_size) { // Used to allocate memory for dynamic arrays
    return malloc(target_size + amnt * child_size);
}

#define RED 0xF800;
#define GREY 0x8080;
const int SW = 1920;
const int SH = 1080;
int SW2 = SW/2;
int SH2 = SH/2;
int SW1 = SW-1;
int SH1 = SH-1;

float sn[360];
float cs[360];

Level *level = NULL;

/*inline void drawWall(Uint32 *pixels, float sx0, float sx1, float sy0, float sy1, float sy2, float sy3) {
    float dx = sx1 - sx0;
    if (dx < 1) {
        dx = 1;
    }
    for (int x=sx0;x<sx1;x++) {
        if (x > 1 & x < SW1) {
            float t = (x-sx0)/dx;
            int Y1 = (int)((1-t) * sy0 + t * sy1);
            int Y2 = (int)((1-t) * sy2 + t * sy3);

            if (Y1 < 1) {
                Y1 = 1;
            } else if (Y1 > SH1) {
                Y1 = SH1;
            }
            if (Y2 < 1) {
                Y2 = 1;
            } else if (Y2 > SH1) {
                Y2 = SH1;
            }
            int index = Y1*SW+x;
            for (int y=Y1;y<Y2;y++) {
                pixels[index] = RED;
                index += SW;
            }
        }
    }
}*/

//void drawPixel(Uint32 *pixels, int x, int y) {
//    if (x > 1 & x < SW1 & y > 1 & y < SH1) {
//        pixels[y*SW+x] = RED;
//    }
//}
int clamp(int value, int minVal, int maxVal) {
    return (value < minVal) ? minVal : (value > maxVal) ? maxVal : value;
}

void drawWall(Uint16 *pixels, float sx0, float sx1, 
    float sy0, float sy1, float sy2, float sy3, 
    portalCull portalBounds, int *ceilingLut, int *floorLut, int flat) {

float dx = (sx1 > sx0) ? (sx1 - sx0) : 1.0f; 
int startX = clamp((int)sx0, portalBounds.x0, portalBounds.x1);
int endX   = clamp((int)sx1, portalBounds.x0, portalBounds.x1);
startX = clamp(startX, 1, SW1);
endX   = clamp(endX, 1, SW1);

int numX = endX - startX;
if (numX <= 0)
return;  // Early exit if nothing to draw.

int topBoundary[numX];
int bottomBoundary[numX];

float invDx = 1.0f / dx;
float invCullDx = 1.0f / portalBounds.dxCull;

for (int x = startX; x < endX; x++) {
int idx = x - startX;
float t = (x - sx0) * invDx;
float st = (x - portalBounds.x0) * invCullDx;

int yTop    = (int)((1.0f - t) * sy0 + t * sy1);
int yBottom = (int)((1.0f - t) * sy2 + t * sy3);
if (yTop > yBottom) { 
  int temp = yTop; 
  yTop = yBottom; 
  yBottom = temp; 
}

int clampY0 = (int)((1.0f - st) * portalBounds.y0 + st * portalBounds.y1);
int clampY1 = (int)((1.0f - st) * portalBounds.y2 + st * portalBounds.y3);
yTop    = clamp(yTop, clampY0, clampY1);
yBottom = clamp(yBottom, clampY0, clampY1);
yTop    = clamp(yTop, 1, SH1);
yBottom = clamp(yBottom, 1, SH1);

topBoundary[idx]    = yTop;
bottomBoundary[idx] = yBottom;

if (flat == 0) {
  ceilingLut[x] = yTop;
  floorLut[x]   = yBottom;
} else if (flat == 1) {
  ceilingLut[x] = yTop;
} else if (flat == 2) {
  floorLut[x] = yBottom;
}
}

// Draw each vertical column using pointer arithmetic
for (int x = startX; x < endX; x++) {
int idx = x - startX;
int yTop = topBoundary[idx];
int yBottom = bottomBoundary[idx];
// Get starting pointer: each column is spaced by SW elements
Uint16 *p = pixels + x + yTop * SW;
int count = yBottom - yTop + 1;
for (int i = 0; i < count; i++) {
  *p = RED;
  p += SW;  // move down one row (SW is the screen width)
}
}
}

void draw_flat(Uint16 *pixels, int *lut, int flat, portalCull portalBounds) {
    int X0 = clamp(clamp((int)0, portalBounds.x0, portalBounds.x1), 1, SW1);
    int X1 = clamp(clamp((int)SW1, portalBounds.x0, portalBounds.x1), 1, SW1);

    if (flat == 1) {
        for (int x=X0; x<X1; x++) {
            float st = (x-portalBounds.x0)/portalBounds.dxCull;

            int clampY0 = (int)((1.0f - st) * portalBounds.y0 + st * portalBounds.y1);
            int clampY1 = (int)((1.0f - st) * portalBounds.y2 + st * portalBounds.y3);

            int Y1 = clamp(clamp(lut[x], clampY0, clampY1), 1, SH1);
            int Y0 = clamp(clamp(1, clampY0, clampY1), 1, SH1);
            int row = Y0*SW+x;
            for (int y=Y0; y<Y1; y++) {
                pixels[row] = GREY;
                row = row + SW;
            }
        }
    } else if (flat == 2) {
        for (int x=X0; x<X1; x++) {
            float st = (x-portalBounds.x0)/portalBounds.dxCull;

            int clampY0 = (int)((1.0f - st) * portalBounds.y0 + st * portalBounds.y1);
            int clampY1 = (int)((1.0f - st) * portalBounds.y2 + st * portalBounds.y3);

            int Y1 = clamp(clamp(SH1, clampY0, clampY1), 1, SH1);
            int Y0 = clamp(clamp(lut[x], clampY0, clampY1), 1, SH1);
            int row = Y0*SW+x;
            for (int y=Y0; y<Y1; y++) {
                pixels[row] = GREY;
                row = row + SW;
            }
        }
    }
}
inline void clip(float *ax, float *ay, float bx, float by, float px1, float py1, float px2, float py2) {
    float a = (px1 - px2) * (*ay - py2) - (py1 - py2) * (*ax - px2);
    float b = (py1 - py2) * (*ax - bx) - (px1 - px2) * (*ay - by);
    float t = a / (b+1);
    
    *ax = *ax - t * (bx - *ax);
    *ay = *ay - t * (by - *ay);
}

int pointWallCollision(float px, float py, float ax, float ay, float bx, float by, float radius) {
    float abx = bx-ax;
    float aby = by-ay;
    float apx = px-ax;
    float apy = py-ay;

    float ab2 = abx*abx + aby*aby;
    float ap_ab = apx*abx + apy*aby;
    float t = ap_ab/ab2;

    float cX = ax+t*abx;
    float cY = ay+t*aby;

    float dx = px-cX;
    float dy = py-cY;
    float colliding = (dx*dx + dy*dy);
    if (colliding < (radius*radius)) {
        return 1;
    } else {
        return 0;
    }
}
void drawSector(Uint16 *pixels, sector *Sector, player character, float pSn, float pCs, portalCull portalBounds) {
    wall *walls = Sector->walls;
    size_t wallCount = Sector->count;
    
    float px = character.x;
    float py = character.y;
    float pz = character.z;
    float wz0 = pz - Sector->height - Sector->elevation;
    float wz1 = pz - Sector->elevation;

    int ceilingLut[SW1];
    int floorLut[SW1];
    
    portalRender *portalQueue = NULL;
    int portalCapacity = 0;
    int portalCount = 0;

    for (int wallIndex=0;wallIndex<wallCount;wallIndex++) {
        wall cWall = walls[wallIndex];
        float rx0 = cWall.x0 - px;
        float ry0 = cWall.y0 - py;
        float rx1 = cWall.x1 - px;
        float ry1 = cWall.y1 - py;

        float tx0 = rx0 * pCs - ry0 * pSn;
        float ty0 = ry0 * pCs + rx0 * pSn;
        float tx1 = rx1 * pCs - ry1 * pSn;
        float ty1 = ry1 * pCs + rx1 * pSn;
        if (ty0 < 1 & ty1 < 1 & !cWall.is_portal) {
            continue;
        }
        if (((((cWall.y1-cWall.y0) * rx1) + (-(cWall.x1-cWall.x0)*ry1)) > -1)) {
            continue;
        }
        int hasClipped = 0;
        if (ty0 < 1) {
            clip(&tx0, &ty0, tx1, ty1, 1.0, 1.0, (float)SW, 1.0);
            hasClipped = 1;
        } 
        if (ty1 < 1) {
            clip(&tx1, &ty1, tx0, ty0, 1.0, 1.0, (float)SW, 1.0);
            hasClipped = 2;
        }
        float inv_ty0 = (1/ty0)*character.focalLength;
        float inv_ty1 = (1/ty1)*character.focalLength;

        if (cWall.is_portal == 0) {
            float sx0 = tx0 * inv_ty0 + SW2;
            float sx1 = tx1 * inv_ty1 + SW2;

            float sy0 = wz0 * inv_ty0 + SH2;
            float sy1 = wz0 * inv_ty1 + SH2;
            float sy2 = wz1 * inv_ty0 + SH2;
            float sy3 = wz1 * inv_ty1 + SH2;
            // We can finally draw our wall!
            //drawWallASM(pixels, sx0, sx1, sy0, sy1, sy2, sy3);
            drawWall(pixels, sx0, sx1, sy0, sy1, sy2, sy3, portalBounds, ceilingLut, floorLut, 0);
        } else {
            float pz1 = wz0 + Sector->height;
            float pz0 = pz1 - cWall.portal_bottom;
            float pz2 = wz0;
            float pz3 = pz2 + cWall.portal_top;

            // Perspective
            float sx0 = tx0 * inv_ty0 + SW2;
            float sy0 = pz0 * inv_ty0 + SH2;

            float sx1 = tx1 * inv_ty1 + SW2;
            float sy1 = pz0 * inv_ty1 + SH2;

            float sy2 = pz1 * inv_ty0 + SH2;
            float sy3 = pz1 * inv_ty1 + SH2;

            float sy4 = pz2 * inv_ty0 + SH2;
            float sy5 = pz2 * inv_ty1 + SH2;

            float sy6 = pz3 * inv_ty0 + SH2;
            float sy7 = pz3 * inv_ty1 + SH2;
            
            drawWall(pixels, sx0, sx1, sy0, sy1, sy2, sy3, portalBounds, ceilingLut, floorLut, 2);
            drawWall(pixels, sx0, sx1, sy4, sy5, sy6, sy7, portalBounds, ceilingLut, floorLut, 1);

            if (portalCount >= portalCapacity) {
                int newCapacity = (portalCapacity == 0) ? 1 : portalCapacity * 2;
                portalRender *temp = realloc(portalQueue, newCapacity * sizeof(portalRender));
                if (!temp) {
                    perror("Failed to realloc memory for portal queue");
                    free(portalQueue);
                    exit(EXIT_FAILURE);
                }
                portalQueue = temp;
                portalCapacity = newCapacity;
            }

            portalRender newPortal;
            newPortal.clipped = 0;
            newPortal.sector_link = cWall.portal_link;
            
            if (hasClipped > 0) {
                int portal_within_screen = 0;
                if ((sx0 >= 0 & sx0 <= SW) || (sx1 >= 0 & sx1 <= SW) || (sx0 < 0 & sx1 > SW) || (sx1 < 0 & sx0 > SW)) {
                    portal_within_screen = 1;
                }
                if (portal_within_screen == 1) {
                    int collision = pointWallCollision(character.x, character.y, cWall.x0, cWall.y0, cWall.x1, cWall.y1, 2.0);
                    if (collision == 1) {
                        if (hasClipped == 1) {
                            sx0 = 1;
                        } else if (hasClipped == 2) {
                            sx1 = SW;
                        }
                    }
                } else {
                    int collision = pointWallCollision(character.x, character.y, cWall.x0, cWall.y0, cWall.x1, cWall.y1, 2.0);
                    if (collision == 1) {
                        sx0 = 1;
                        sx1 = SW;
                    }
                }
            }
            portalCull newPortalBounds = {sx0, sx1, sy6, sy7, sy0, sy1, sx1-sx0};
            newPortal.portalBounds = newPortalBounds;
            portalQueue[portalCount++] = newPortal;

        }
    }
    draw_flat(pixels, ceilingLut, 1, portalBounds);
    draw_flat(pixels, floorLut, 2, portalBounds);
    for (int portalIndex=0; portalIndex<portalCapacity; portalIndex++) {
        portalRender portalInfo = portalQueue[portalIndex];
        sector *portalSector = level->sectors[portalInfo.sector_link];
        if (portalInfo.sector_link < 0 || portalInfo.sector_link > level->count) {
            perror("Failed to locate sector for portal traversal");
            exit(EXIT_FAILURE);
        }
        drawSector(pixels, portalSector, character, pSn, pCs, portalInfo.portalBounds);
    }
    free(portalQueue);
}

void beginDrawSector(Uint16 *pixels, sector *Sector, player character, float pSn, float pCs, portalCull portalBounds) {
    drawSector(pixels, Sector, character, pSn, pCs, portalBounds);
}

void testDraw(Uint16 *pixels, int pitch) {
    for (int y = 0; y < SH; y++) {
        for (int x = 0; x < SW; x++) {
            pixels[y * pitch + x] = RED;
        }
    }
}

void clearPixels(Uint16 *pixels, int pitch) {
    for (int y=0; y <SH; y++) {
        for (int x=0; x<SW; x++) {
            pixels[y*pitch+x] = 0;
        }
    }
}

int avg = 0;
int avgI = 0;

int isPlayerWithinConvexSector(player character) {
    int i = 0;
    
    for (int s=0; s<level->count; s++) {
        sector *Sector = level->sectors[s];
        
        int inside = 1;
        wall *walls = Sector->walls;
        size_t wallCount = Sector->count;
        
        for (int w=0; w<wallCount; w++) {
            wall cWall = walls[w];
            float cross = (cWall.x1-cWall.x0) * (character.y-cWall.y0) - (cWall.y1-cWall.y0) * (character.x-cWall.x0);
            if (cross > 0) {
                inside = 0;
            }
        }
        
        if (inside == 1) {
            return i;
        }
        i ++;
    }
    return -1;
}

// seems to be that clipping or something related to rendering other sectors is harming performance

int main(int argc, char* argv[]) {
    // Let's create our level
    level = malloc(sizeof(Level) + 1 * sizeof(sector *));
    if (!level) {
        perror("Failed to allocate level");
    }
    level->count = 2;
    //
    level->sectors[0] = malloc(sizeof(sector) + 4 * sizeof(wall));
    sector *sect = level->sectors[0];
    sect->count = 4;
    sect->height = 10.0;
    sect->elevation = 0.0;
    sect->walls[0] = (wall){40, 0, 0, 0, 1, 2, 2, 1};
    sect->walls[1] = (wall){0, 0, 0, 40, 0, 0, 0, 0};
    sect->walls[2] = (wall){40, 40, 40, 0, 0, 0, 0, 0};
    sect->walls[3] = (wall){0, 40, 40, 40, 0, 0, 0, 0};
    //
    level->sectors[1] = malloc(sizeof(sector) + 4 * sizeof(wall));
    sector *sect1 = level->sectors[1];
    sect1->count = 4;
    sect1->height = 6.0;
    sect1->elevation = 2.0;
    sect1->walls[0] = (wall){40, 0, -5, -5, 0, 0, 0, 0};
    sect1->walls[1] = (wall){0, 0, -5, 0, 0, 0, 0, 0};
    sect1->walls[2] = (wall){40, 40, 0, -5, 0, 0, 0, 0};
    sect1->walls[3] = (wall){0, 40, 0, 0, 0, 0, 0, 0};

    // Let's setup our sine and cosine lookup tables
    for (int i = 0; i < 360; i++) {
        double a = (double)(i)*M_PI/180;
        sn[i] = sin(a);
        cs[i] = cos(a);
    }
    SDL_Window *window = SDL_CreateWindow("SDL2 Software Renderer Example",
                                          SDL_WINDOWPOS_CENTERED,
                                          SDL_WINDOWPOS_CENTERED,
                                          SW, SH, 0);
    if (!window) {
        fprintf(stderr, "Window creation failed: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    if (!renderer) {
        fprintf(stderr, "Renderer creation failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    SDL_Surface *surface = SDL_GetWindowSurface(window);
    SDL_Surface *rgb565_surface = SDL_ConvertSurfaceFormat(surface, SDL_PIXELFORMAT_RGB565, 0);
    Uint16 *pixels = (Uint16 *)rgb565_surface->pixels;
    //int pitch = surface->pitch / sizeof(Uint32);
    // Main loop: wait for the user to close the window
    SDL_Event event;
    int quit = 0;

    //sector *currentSector = level->sectors[0];
    player character = {20.0, 20.0, 5.0, 0.0, 180.0, 0.0, 90.0};
    // Let's calculate our focalLength value
    float f = character.fov*M_PI/180;
    float tanFOV = tan(f/2);
    float focalLength = SW2/tanFOV;
    character.focalLength = focalLength;

    portalCull portalBounds = {1.0, (float)SW1, 1.0, 1.0, SH1, SH1, SW1-1.0};

    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            } else if (event.type == SDL_KEYDOWN) {
                float a = (character.yaw)*M_PI/180;
                if (event.key.keysym.sym == SDLK_w) {
                    character.x = character.x + sin(a);
                    character.y = character.y + cos(a);
                } else if (event.key.keysym.sym == SDLK_s) {
                    character.x = character.x - sin(a);
                    character.y = character.y - cos(a);
                }
                if (event.key.keysym.sym == SDLK_a) {
                    character.yaw -= 1;
                } else if (event.key.keysym.sym == SDLK_d) {
                    character.yaw += 1;
                }
            }
        }
        // Get the performance frequency (ticks per second)
        double frequency = (double)SDL_GetPerformanceFrequency();

        // Start time
        //SDL_LockSurface(surface);
        //testDraw(pixels, pitch);
        //SDL_UnlockSurface(surface);
        
        // inputs
        if (character.yaw < 0) {
            character.yaw = 360;
        } else if (character.yaw > 360) {
            character.yaw = 0;
        }
        int yaw = ((int)character.yaw) % 360;
        float pSn = sn[yaw];
        float pCs = cs[yaw];
        //SDL_LockSurface(rgb565_surface);
        int sectorIndex = isPlayerWithinConvexSector(character);
        sector *playerSector = level->sectors[sectorIndex];

        Uint64 st = SDL_GetPerformanceCounter(); 
        beginDrawSector(pixels, playerSector, character, pSn, pCs, portalBounds);
        SDL_BlitSurface(rgb565_surface, NULL, surface, NULL);
        SDL_UpdateWindowSurface(window);
        Uint64 et = SDL_GetPerformanceCounter(); 
        //SDL_UnlockSurface(rgb565_surface);
       //memset(pixels, 0, SW*SH*sizeof(Uint16));
        double ft = (et - st) / frequency;
        int fps = (int)(1/ft);
        avg = avg + fps;
        avgI ++;
        if (avgI == 120) {
            char str[10];
            sprintf(str, "%d", avg/avgI);
            SDL_SetWindowTitle(window, str);
            avg = 0;
            avgI = 0;
        }
    }

    // Clean up resources before exiting
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    free(level->sectors[0]);
    free(level);
    return 0;
}