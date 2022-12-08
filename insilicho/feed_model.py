from insilicho import reactor


def F(time, V, bolus_size=0.03, num_bolus=1, bolus_frequency=24):
    def dirac(t):
        # https://stackoverflow.com/questions/63230568/numerical-solution-to-a-differential-equation-containing-a-dirac-delta-function
        def tentfunc(t):
            return max(0, 1 - abs(t))

        N = 10.0
        return sum(
            N * tentfunc(N * (t - (i + 1) * bolus_frequency)) for i in range(num_bolus)
        )

    if 0 < V <= reactor.Vessel.volume:
        return bolus_size * dirac(time)  # L/h
    return 0
