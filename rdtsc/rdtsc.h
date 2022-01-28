#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long x;
    __asm__ volatile ("rdtsc" : "=A" (x));
    return x;
}

#elif defined(__x86_64__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long a, d;
    __asm__ volatile ("rdtsc" : "=a"(a), "=d"(d));
    return (d<<32) | a;
}
#endif
