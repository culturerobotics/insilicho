from insilicho import run


# Objective is final mAbs concentration
def score(exp_res):
    return exp_res["Cmab"][-1]


def run_exp(factor_settings, model_obj=None, Xv=0.3e5, plot=False):
    if not model_obj:
        # Run for 96 hours, sampling every 12 hours
        model_obj = run.GrowCHO({"parameters": {"Ndays": 4, "Nsamples": 2}}, None, None)

    """Run an experiment with the supplied factor settings (falling back to defaults for missing settings)."""

    default_settings = {
        "batch_glc": 100,
        "batch_gln": 100,
        "batch_pH": 7,
        "feed_glc": 300,
        "feed_gln": 300,
        "prod_start_eft": 48,
        "batch_temp": 36,
        "prod_temp": 36,
        "day_0_feed": 0,
        "day_1_feed": 0,
        "day_2_feed": 0,
        "day_3_feed": 0,
    }

    settings = default_settings | factor_settings

    def feed(t, _V):
        feed_array = [
            settings["day_0_feed"],
            settings["day_1_feed"],
            settings["day_2_feed"],
            settings["day_3_feed"],
            settings["day_3_feed"],
        ]  # Repeat last entry to avoid index out of bounds
        return feed_array[int(t // 24)]

    def temp(t):
        if t < settings["prod_start_eft"]:
            return settings["batch_temp"]
        else:
            return settings["prod_temp"]

    model_obj.params.Cglc_feed = settings["feed_glc"]
    model_obj.params.Cgln_feed = settings["feed_gln"]

    model_obj.feed_fn = feed
    model_obj.temp_fn = temp

    res = model_obj.execute(
        {
            "Cglc": settings["batch_glc"],
            "Cgln": settings["batch_gln"],
            "pH": settings["batch_pH"],
            "V": 0.12,
            "Xv": Xv,
        },
        plot=plot,
    )

    return res, score(res)


ranges = {
    "batch_glc": (50, 150),
    "batch_gln": (50, 150),
    "batch_pH": (6.7, 7.4),
    "feed_glc": (100, 500),
    "feed_gln": (100, 500),
    "prod_start_eft": (24, 72),
    "batch_temp": (35, 37),
    "prod_temp": (35, 37),
    "day_0_feed": (0, 0.005),
    "day_1_feed": (0, 0.005),
    "day_2_feed": (0, 0.005),
    "day_3_feed": (0, 0.005),
}
