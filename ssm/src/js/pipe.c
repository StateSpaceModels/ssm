#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#define ParentRead      read_pipe[0]
#define ParentWrite     write_pipe[1]
#define ChildRead       write_pipe[0]
#define ChildWrite      read_pipe[1]

void error(char *msg)
{
    fprintf(stderr, "%s: %s\n", msg, strerror(errno));
    exit(1);
}

int main()
{
    /** Pipe for reading for subprocess */
    int read_pipe[2];
    if (pipe(read_pipe) == -1) {
	error("Can't create the pipe");
    }

    /** Pipe for writing to subprocess */
    int write_pipe[2];
    if (pipe(write_pipe) == -1) {
	error("Can't create the pipe");
    }

    pid_t pid = fork();
    if (pid == -1) {
	error("Can't fork process");
    } 

    if (!pid) { //Child process

	close(ParentRead);
	close(ParentWrite);
	dup2(ChildRead, 0);
	dup2(ChildWrite, 1);


	if (execlp("node", "node", "./pipe.js", NULL) == -1) {
	    error("Can't run script");
	}
    }
    
    close(ChildRead);
    close(ChildWrite);

    dup2(ParentRead, 0);
    
    printf("starting...\n");

    int i;
    for(i=0; i<2;i++){
    
	write(ParentWrite, "abc", 3);

	char line[255];    
	if (fgets(line, 255, stdin) != NULL) {
	    printf("received fgets: %s", line);
	}

	if (fgets(line, 255, stdin) != NULL) {
	    printf("received fgets: %s", line);
	}

    }

    return 0;
}
