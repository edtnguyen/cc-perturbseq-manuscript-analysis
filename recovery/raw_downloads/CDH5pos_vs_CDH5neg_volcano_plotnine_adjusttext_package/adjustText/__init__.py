"""Deterministic adjustText-compatible automatic label placement."""
from __future__ import annotations

from time import perf_counter
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch

__version__ = "1.4.0-compatible-vendored"


def _pair(value, default=(0.0, 0.0)):
    if value is None:
        return float(default[0]), float(default[1])
    if np.isscalar(value):
        return float(value), float(value)
    values = tuple(value)
    if len(values) != 2:
        raise ValueError(f"Expected scalar or length-2 value, got {value!r}")
    return float(values[0]), float(values[1])


def _text_bboxes(texts, renderer, expand=(1.0, 1.0)):
    return np.asarray(
        [
            (
                box.xmin,
                box.xmax,
                box.ymin,
                box.ymax,
            )
            for box in (
                text.get_window_extent(renderer=renderer).expanded(*expand)
                for text in texts
            )
        ],
        dtype=float,
    )


def _centers(boxes):
    return np.column_stack(
        ((boxes[:, 0] + boxes[:, 1]) / 2, (boxes[:, 2] + boxes[:, 3]) / 2)
    )


def _movement_mask(spec, axis, values):
    if axis not in spec:
        return np.zeros_like(values)
    out = values.copy()
    if f"{axis}+" in spec:
        out[out < 0] = 0
    elif f"{axis}-" in spec:
        out[out > 0] = 0
    return out


def _set_display_centers(texts, centers, transform):
    positions = transform.inverted().transform(centers)
    for text, position in zip(texts, positions):
        text.set_horizontalalignment("center")
        text.set_verticalalignment("center")
        text.set_position(position)


def _orientation(p, q, r):
    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])


def _segments_cross(a, b):
    p1, p2 = a[:2], a[2:]
    q1, q2 = b[:2], b[2:]
    o1 = _orientation(p1, p2, q1)
    o2 = _orientation(p1, p2, q2)
    o3 = _orientation(q1, q2, p1)
    o4 = _orientation(q1, q2, p2)
    return (o1 * o2 < -1e-9) and (o3 * o4 < -1e-9)


def adjust_text(
    texts,
    x=None,
    y=None,
    objects=None,
    target_x=None,
    target_y=None,
    avoid_self=True,
    prevent_crossings=True,
    force_text=(0.1, 0.2),
    force_static=(0.1, 0.2),
    force_pull=(0.01, 0.01),
    force_explode=(0.1, 0.5),
    pull_threshold=10,
    expand=(1.05, 1.2),
    max_move=(10, 10),
    explode_radius="auto",
    ensure_inside_axes=True,
    expand_axes=False,
    only_move=None,
    ax=None,
    min_arrow_len=5,
    time_lim=None,
    iter_lim=None,
    *args,
    **kwargs,
):
    if not texts:
        return [], []

    ax = ax or plt.gca()
    fig = ax.figure
    try:
        fig.draw_without_rendering()
    except AttributeError:
        fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    transform = texts[0].get_transform()

    initial_data = np.asarray(
        [text.get_unitless_position() for text in texts], dtype=float
    )
    if target_x is not None and target_y is not None:
        target_data = np.column_stack(
            (np.asarray(target_x, float), np.asarray(target_y, float))
        )
    else:
        target_data = initial_data.copy()
    if len(target_data) != len(texts):
        raise ValueError("target_x/target_y must have the same length as texts")
    target_display = transform.transform(target_data)

    point_display = np.empty((0, 2), dtype=float)
    if x is not None and y is not None:
        point_display = transform.transform(
            np.column_stack((np.asarray(x, float), np.asarray(y, float)))
        )
    if avoid_self:
        point_display = np.vstack((point_display, target_display))

    force_text = _pair(force_text, (0.1, 0.2))
    force_static = _pair(force_static, (0.1, 0.2))
    force_pull = _pair(force_pull, (0.01, 0.01))
    force_explode = _pair(force_explode, (0.1, 0.5))
    max_move = _pair(max_move, (10, 10)) if max_move is not None else (np.inf, np.inf)
    expand = _pair(expand, (1.05, 1.2))

    only_move = only_move or {
        "text": "xy",
        "static": "xy",
        "explode": "xy",
        "pull": "xy",
    }
    if isinstance(only_move, str):
        only_move = {key: only_move for key in ("text", "static", "explode", "pull")}
    for key in ("text", "static", "explode", "pull"):
        only_move.setdefault(key, "xy")

    boxes = _text_bboxes(texts, renderer)
    centers = _centers(boxes)
    widths = boxes[:, 1] - boxes[:, 0]
    heights = boxes[:, 3] - boxes[:, 2]

    if explode_radius == "auto":
        explode_radius = max(float(np.mean(widths)), float(np.mean(heights)))

    rng = np.random.default_rng(42)
    jitter = rng.uniform(-1.0, 1.0, size=centers.shape)
    jitter[:, 0] *= float(explode_radius) * force_explode[0]
    jitter[:, 1] *= float(explode_radius) * force_explode[1]
    jitter[:, 0] = _movement_mask(only_move["explode"], "x", jitter[:, 0])
    jitter[:, 1] = _movement_mask(only_move["explode"], "y", jitter[:, 1])
    centers += jitter
    _set_display_centers(texts, centers, transform)

    if iter_lim is None:
        iter_lim = 600
    start = perf_counter()
    axes_box = ax.patch.get_window_extent(renderer=renderer)
    stable_rounds = 0

    for _ in range(int(iter_lim)):
        if time_lim is not None and perf_counter() - start >= float(time_lim):
            break

        boxes = _text_bboxes(texts, renderer, expand=expand)
        centers = _centers(boxes)
        movement = np.zeros_like(centers)
        overlap_count = 0

        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                overlap_x = min(boxes[i, 1], boxes[j, 1]) - max(boxes[i, 0], boxes[j, 0])
                overlap_y = min(boxes[i, 3], boxes[j, 3]) - max(boxes[i, 2], boxes[j, 2])
                if overlap_x <= 0 or overlap_y <= 0:
                    continue

                overlap_count += 1
                dx = centers[i, 0] - centers[j, 0]
                dy = centers[i, 1] - centers[j, 1]

                if overlap_x < overlap_y:
                    sign = np.sign(dx) or (1.0 if (i + j) % 2 == 0 else -1.0)
                    amount = (overlap_x / 2 + 1.25) * max(force_text[0] * 6.0, 0.65)
                    movement[i, 0] += sign * amount
                    movement[j, 0] -= sign * amount
                else:
                    sign = np.sign(dy) or (1.0 if (i + j) % 2 == 0 else -1.0)
                    amount = (overlap_y / 2 + 1.25) * max(force_text[1] * 4.0, 0.65)
                    movement[i, 1] += sign * amount
                    movement[j, 1] -= sign * amount

        static_hits = 0
        if len(point_display):
            for i, box in enumerate(boxes):
                inside = (
                    (point_display[:, 0] >= box[0])
                    & (point_display[:, 0] <= box[1])
                    & (point_display[:, 1] >= box[2])
                    & (point_display[:, 1] <= box[3])
                )
                for point in point_display[inside]:
                    static_hits += 1
                    dx = centers[i, 0] - point[0]
                    dy = centers[i, 1] - point[1]
                    if abs(dy) <= abs(dx) * 1.4:
                        sign = np.sign(dy) or 1.0
                        movement[i, 1] += sign * max(force_static[1] * 18.0, 1.0)
                    else:
                        sign = np.sign(dx) or 1.0
                        movement[i, 0] += sign * max(force_static[0] * 18.0, 1.0)

        pull = target_display - centers
        pull_x = np.where(
            np.abs(pull[:, 0]) > pull_threshold,
            pull[:, 0] * force_pull[0],
            0.0,
        )
        pull_y = np.where(
            np.abs(pull[:, 1]) > pull_threshold,
            pull[:, 1] * force_pull[1],
            0.0,
        )
        movement[:, 0] += _movement_mask(only_move["pull"], "x", pull_x)
        movement[:, 1] += _movement_mask(only_move["pull"], "y", pull_y)

        movement[:, 0] = _movement_mask(only_move["text"], "x", movement[:, 0])
        movement[:, 1] = _movement_mask(only_move["text"], "y", movement[:, 1])
        movement[:, 0] = np.clip(movement[:, 0], -max_move[0], max_move[0])
        movement[:, 1] = np.clip(movement[:, 1], -max_move[1], max_move[1])

        new_centers = centers + movement
        if ensure_inside_axes:
            half_w = (boxes[:, 1] - boxes[:, 0]) / 2
            half_h = (boxes[:, 3] - boxes[:, 2]) / 2
            new_centers[:, 0] = np.clip(
                new_centers[:, 0],
                axes_box.xmin + half_w + 1,
                axes_box.xmax - half_w - 1,
            )
            new_centers[:, 1] = np.clip(
                new_centers[:, 1],
                axes_box.ymin + half_h + 1,
                axes_box.ymax - half_h - 1,
            )

        _set_display_centers(texts, new_centers, transform)
        max_step = float(np.max(np.abs(movement))) if movement.size else 0.0
        if overlap_count == 0 and static_hits == 0 and max_step < 0.35:
            stable_rounds += 1
            if stable_rounds >= 4:
                break
        else:
            stable_rounds = 0

    if prevent_crossings and len(texts) > 1:
        for _ in range(5):
            boxes = _text_bboxes(texts, renderer)
            centers = _centers(boxes)
            changed = False
            for i in range(len(texts)):
                for j in range(i + 1, len(texts)):
                    seg_i = np.r_[centers[i], target_display[i]]
                    seg_j = np.r_[centers[j], target_display[j]]
                    if not _segments_cross(seg_i, seg_j):
                        continue
                    current = np.linalg.norm(centers[i] - target_display[i]) + np.linalg.norm(
                        centers[j] - target_display[j]
                    )
                    swapped = np.linalg.norm(centers[j] - target_display[i]) + np.linalg.norm(
                        centers[i] - target_display[j]
                    )
                    if swapped <= current * 1.10:
                        centers[[i, j]] = centers[[j, i]]
                        _set_display_centers(texts, centers, transform)
                        changed = True
            if not changed:
                break

    arrowprops = dict(kwargs.pop("arrowprops", {}) or {})
    patches = []
    final_centers = _centers(_text_bboxes(texts, renderer))

    for text, center, target_data_i, target_display_i in zip(
        texts, final_centers, target_data, target_display
    ):
        if not arrowprops or np.linalg.norm(center - target_display_i) < min_arrow_len:
            continue
        props = dict(arrowprops)
        props.setdefault("shrinkA", 2.0)
        props.setdefault("shrinkB", 4.0)
        props.setdefault("connectionstyle", "arc3,rad=0")
        patch = FancyArrowPatch(
            posA=text.get_unitless_position(),
            posB=target_data_i,
            patchA=text,
            transform=transform,
            *args,
            **kwargs,
            **props,
        )
        ax.add_patch(patch)
        patches.append(patch)

    return list(texts), patches
