"""Shared test doubles and design helpers."""
from __future__ import annotations

from typing import NamedTuple

from speaker.drivers import Driver, load_driver
from speaker.layout import Layout, compute_layout
from speaker.params import Params, load_params


class DesignTuple(NamedTuple):
    p: Params
    coax: Driver
    bass: Driver

    @property
    def layout(self) -> Layout:
        return compute_layout(self.p, self.coax, self.bass)


def default_design() -> DesignTuple:
    """params/default.toml with the drivers it names."""
    p = load_params()
    return DesignTuple(p, load_driver(p.coax), load_driver(p.bass))
