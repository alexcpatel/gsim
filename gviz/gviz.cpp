#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>

#define DIST_SCALE 20.0e15 // meters
#define INIT_MASS 4.0e30 // kg

typedef struct {
  double m;
  double x;
  double y;
} pos_t;

static void error_callback(int error, const char* description) {
  fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key,
                         int scancode, int action, int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
      glfwSetWindowShouldClose(window, GL_TRUE);
}

static int setup(char *pos_path, pos_t ***positions,
                 size_t *num_bodies, size_t *num_steps) {
  /* get header line */
  std::ifstream infile(pos_path);
  if (!(infile >> *num_bodies >> *num_steps)) return -1;
  (*num_steps)++;

  /* allocate memory */
  size_t si, bi;
  *positions = (pos_t **)malloc(*num_steps * sizeof(pos_t *));
  if (*positions == NULL) return -1;
  for (si = 0; si < *num_steps; si++) {
    (*positions)[si] = (pos_t *)malloc(*num_bodies * sizeof(pos_t));
    if ((*positions)[si] == NULL) return -1;
  }

  /* read file and fill memory */
  const double dist_divisor = 2.0 / DIST_SCALE;
  const double mass_divisor = 6.0 / INIT_MASS;
  double m, x, y;
  pos_t pos;
  while (infile >> si >> bi >> m >> x >> y) {
    pos.m = m * mass_divisor;
    pos.x = x * dist_divisor - 1.0;
    pos.y = y * dist_divisor - 1.0;
    ((*positions)[si])[bi] = pos;
  }

  return 0;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: ./gviz <path to position file>\n");
    exit(EXIT_FAILURE);
  }

  pos_t **positions;
  size_t num_bodies, num_steps;
  if (setup(argv[1], &positions, &num_bodies, &num_steps)) {
    fprintf(stderr, "Invalid position file.\n");
    exit(EXIT_FAILURE);
  }

  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) exit(EXIT_FAILURE);
  window = glfwCreateWindow(640, 480, "gsim visualizer", NULL, NULL);
  if (!window) {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);

  size_t si = 0, bi;
  pos_t *step_positions, pos;

  while (!glfwWindowShouldClose(window)) {
    /* initialize window and position space */
    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float) height;

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
    glMatrixMode(GL_MODELVIEW);

    /* render body positions */
    step_positions = positions[si];
    for (bi = 0; bi < num_bodies; bi++) {
      pos = step_positions[bi];

      /* set point options */
      glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
      glEnable(GL_POINT_SMOOTH);
      glLoadIdentity();
      glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
      glPointSize(pos.m);
      glBegin(GL_POINTS);
      glColor3f(1.f, 1.f, 1.f);

      /* draw vertex */
      glVertex2f(pos.x, pos.y);

      glEnd();
    }

    /* end opengl rendering */
    
    glfwSwapBuffers(window);
    glfwPollEvents();

    /* increment step and continue */
    si++;
    if (si >= num_steps) si = 0;
    usleep(1000); // sleep 100 ms
  }

  glfwDestroyWindow(window);
  glfwTerminate();

  /* free memory */
  for (si = 0; si < num_steps; si++) free(positions[si]);
  free(positions);

  exit(EXIT_SUCCESS);
}