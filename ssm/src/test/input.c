#include <stdio.h>
#include <jansson.h>

int main(int argc, char *argv[])
{
    json_error_t error;
    json_t *data = json_loadf(stdin, 0, &error);
    if(!data) {
        printf("error: %s\n", error.text);
        exit(EXIT_FAILURE);
    }
    printf("read: %s", json_dumps(data, 0));

    return 0;
}
