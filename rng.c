#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>

double get_random_double(void) {
    uint32_t r;

#if defined(_WIN32) || defined(_WIN64)
    // Windows: Use modern BCryptGenRandom
    #include <windows.h>
    #include <bcrypt.h>
    #pragma comment(lib, "bcrypt.lib")
    
    NTSTATUS status = BCryptGenRandom(NULL, (PUCHAR)&r, sizeof(r), BCRYPT_USE_SYSTEM_PREFERRED_RNG);
    
    if (!BCRYPT_SUCCESS(status)) {
        fprintf(stderr, "BCryptGenRandom failed with error code: 0x%08x\n", status);
        exit(EXIT_FAILURE);
    }

#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    // macOS/BSD: Use arc4random_buf for better efficiency
    arc4random_buf(&r, sizeof(r));

#else // Linux and other Unix-like systems
    // Linux/Unix: Use getrandom() syscall with fallback to /dev/urandom
    #include <fcntl.h>
    #include <unistd.h>
    
    #if defined(__linux__) && defined(SYS_getrandom)
    #include <sys/syscall.h>
    #include <linux/random.h>
    
    // Try getrandom() first (Linux 3.17+)
    if (syscall(SYS_getrandom, &r, sizeof(r), 0) == sizeof(r)) {
        // Successfully got random bytes
    } else {
        // Fallback to /dev/urandom
        int fd = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
        if (fd < 0 || read(fd, &r, sizeof(r)) != sizeof(r)) {
            perror("Failed to get random bytes");
            if (fd >= 0) close(fd);
            exit(EXIT_FAILURE);
        }
        close(fd);
    }
    
    #else
    // Fallback to /dev/urandom on systems without getrandom
    int fd = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    if (fd < 0 || read(fd, &r, sizeof(r)) != sizeof(r)) {
        perror("Failed to get random bytes from /dev/urandom");
        if (fd >= 0) close(fd);
        exit(EXIT_FAILURE);
    }
    close(fd);
    #endif

#endif

    return (double)r / UINT32_MAX;
}

