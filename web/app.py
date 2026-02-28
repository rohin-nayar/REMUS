import sys
import subprocess
import threading
from pathlib import Path
from flask import Flask, render_template, request, jsonify

app = Flask(__name__)

BASE_DIR = Path(__file__).parent
SIM_CPP = BASE_DIR.parent / "simulation" / "REMUS_Code.cpp"
BIN_EXT = ".exe" if sys.platform == "win32" else ""
SIM_BIN = BASE_DIR / f"remus_sim{BIN_EXT}"

sim_lock = threading.Lock()


def compile_simulation():
    """Compile the C++ binary for reference/verification."""
    result = subprocess.run(
        ["g++", "-include", "cmath", str(SIM_CPP), "-o", str(SIM_BIN)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"WARNING: C++ compilation failed:\n{result.stderr}")
        return False
    print(f"C++ binary compiled: {SIM_BIN}")
    return True


# ── Python RK4 simulation ────────────────────────────────────────────────────

def run_simulation_python(env, stages):
    """
    RK4 simulation matching the physics in REMUS_Code.cpp.
    Fixes the v=0 singularity: when v=0 thrust is applied upward (sign = +1).
    Returns list of {t, h, v, m} dicts.
    """
    rho = env["rho"]
    g   = env["g"]
    Cd  = env["Cd"]
    v   = env["v0"]
    h   = env["h0"]

    nstage = len(stages)

    def sign(x):
        return 1.0 if x >= 0 else -1.0

    def burn_time(s):
        return (s["m0"] - s["mR"]) / s["mf"]

    def total_mass(stage_idx):
        """Initial mass for a stage: sum of this stage's m0 + all subsequent stages' m0."""
        return sum(stages[j]["m0"] for j in range(stage_idx, nstage))

    rows = []
    duration = 0.0
    maxh = h

    # ── Powered phase (each stage) ──────────────────────────────────────────
    for q, stage in enumerate(stages):
        m    = total_mass(q)
        mf   = stage["mf"]
        ue   = stage["ue"]
        A    = stage["A"]
        dt   = stage["t"]
        tburn = burn_time(stage)

        drag_cf = 0.5 * rho * Cd * A

        t = 0.0
        while t <= tburn:
            duration += dt

            def vdot(vel, mass):
                drag  = drag_cf * vel * abs(vel) / mass
                thrust = mf * ue * sign(vel) / mass  # sign(v)=+1 on ascent; avoids 0/0
                return -g - drag + thrust

            def hdot(vel):
                return vel

            vk1 = dt * vdot(v, m)
            hk1 = dt * hdot(v)

            vk2 = dt * vdot(v + 0.5 * vk1, m)
            hk2 = dt * hdot(v + 0.5 * vk1)

            vk3 = dt * vdot(v + 0.5 * vk2, m)
            hk3 = dt * hdot(v + 0.5 * vk2)

            vk4 = dt * vdot(v + vk3, m)
            hk4 = dt * hdot(v + vk3)

            m += dt * (-mf)
            v += (vk1 + 2*(vk2 + vk3) + vk4) / 6.0
            h += (hk1 + 2*(hk2 + hk3) + hk4) / 6.0

            if h > maxh:
                maxh = h
            if h < 0:
                h = 0
                rows.append({"t": round(duration, 5), "h": 0.0, "v": round(v, 5), "m": round(m, 5)})
                return rows

            rows.append({"t": round(duration, 5), "h": round(h, 5), "v": round(v, 5), "m": round(m, 5)})
            t += dt

    # ── Freefall phase ──────────────────────────────────────────────────────
    m_final = stages[-1]["mR"]
    A_final = stages[-1]["A"]
    dt      = stages[-1]["t"]
    drag_cf = 0.5 * rho * Cd * A_final

    def vdot_free(vel):
        return -g - drag_cf * vel * abs(vel) / m_final

    while h > 0:
        duration += dt

        vk1 = dt * vdot_free(v)
        hk1 = dt * v

        vk2 = dt * vdot_free(v + 0.5 * vk1)
        hk2 = dt * (v + 0.5 * vk1)

        vk3 = dt * vdot_free(v + 0.5 * vk2)
        hk3 = dt * (v + 0.5 * vk2)

        vk4 = dt * vdot_free(v + vk3)
        hk4 = dt * (v + vk3)

        v += (vk1 + 2*(vk2 + vk3) + vk4) / 6.0
        h += (hk1 + 2*(hk2 + hk3) + hk4) / 6.0

        if h < 0:
            h = 0

        rows.append({"t": round(duration, 5), "h": round(max(h, 0), 5), "v": round(v, 5), "m": round(m_final, 5)})

    return rows


# ── Flask routes ─────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/simulate", methods=["POST"])
def simulate():
    data = request.get_json()
    if not data:
        return jsonify({"error": "No data provided"}), 400

    env    = data.get("env", {})
    stages = data.get("stages", [])

    if not stages:
        return jsonify({"error": "At least one rocket stage is required"}), 400

    try:
        with sim_lock:
            rows = run_simulation_python(env, stages)

        if not rows:
            return jsonify({"error": "Simulation produced no data — check that thrust exceeds gravity"}), 500

        stats = {
            "nStages":         len(stages),
            "maxAltitude":     round(max(r["h"] for r in rows), 2),
            "terminalVelocity": round(rows[-1]["v"], 3),
            "flightTime":      round(rows[-1]["t"], 2),
        }
        return jsonify({"data": rows, "stats": stats})

    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


if __name__ == "__main__":
    compile_simulation()   # compile for reference; simulation runs in Python
    app.run(debug=True, port=5000)
