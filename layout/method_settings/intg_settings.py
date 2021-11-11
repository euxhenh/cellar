import dash_bootstrap_components as dbc
import dash_html_components as html

intg_stvea_settings = dbc.Popover([
    dbc.PopoverHeader("STvEA Settings"),
    dbc.PopoverBody([
        dbc.FormGroup([
            dbc.Label("Clean CODEX"),
            dbc.RadioItems(
                options=[
                    {"label": "True", "value": True},
                    {"label": "False", "value": False},
                ],
                value=True,
                id="intg-clean-codex",
                inline=True
            )
        ]),
        dbc.FormGroup([
            dbc.Label("CODEX clean model"),
            dbc.Select(
                id="intg-clean-codex-method",
                options=[
                    {"label": "gaussian", "value": "gaussian"},
                    {"label": "negative-binomial (slow)", "value": "nb"}
                ],
                value="gaussian"
            )
        ]),
        dbc.FormGroup([
            dbc.Label("Normalize CODEX (used if clean=True and model=nb)"),
            dbc.RadioItems(
                options=[
                    {"label": "True", "value": True},
                    {"label": "False", "value": False},
                ],
                value=False,
                id="intg-clean-codex-normalize",
                inline=True
            )
        ]),
        dbc.FormGroup([
            dbc.Label("Clean CITE protein"),
            dbc.RadioItems(
                options=[
                    {"label": "True", "value": True},
                    {"label": "False", "value": False},
                ],
                value=True,
                id="intg-clean-cite",
                inline=True
            )
        ]),
        dbc.FormGroup([
            dbc.Label("CITE protein clean model"),
            dbc.Select(
                id="intg-clean-cite-method",
                options=[
                    {"label": "gaussian", "value": "gaussian"},
                    {"label": "negative-binomial (slow)", "value": "nb"}
                ],
                value="gaussian"
            )
        ]),
        dbc.FormGroup([
            dbc.Label(
                "Normalize CITE protein (used if clean=True and model=nb)"),
            dbc.RadioItems(
                options=[
                    {"label": "True", "value": True},
                    {"label": "False", "value": False},
                ],
                value=True,
                id="intg-clean-cite-normalize",
                inline=True
            )
        ]),
        html.P([
            "Integrates CITE-seq and CODEX data. Uses an "
            "anchor correction method to map features between "
            "the two. Must be used on a CODEX dataset and "
            "will expand that dataset by adding expresion "
            "levels from the genes that were mapped from "
            "CITE-seq. These genes can be found under "
            "the 'Other' features tab."
        ]),
        html.P([
            "Details: ",
            html.A("https://github.com/CamaraLab/STvEA/"
                   "blob/master/examples/mapping_tutorial.md",
                   href="https://github.com/CamaraLab/"
                   "STvEA/blob/master/examples/mapping_tutorial.md",
                   target="_blank")
        ])
    ])],
    id="intg-STvEA-settings",
    target="intg-STvEA-btn",
    trigger="legacy"
)

intg_stvea_settings_keys = {
    'intg-clean-codex': 'clean_codex',
    'intg-clean-codex-method': 'codex_clean_model',
    'intg-clean-codex-normalize': 'codex_normalize',
    'intg-clean-cite': 'clean_cite',
    'intg-clean-cite-method': 'cite_clean_model',
    'intg-clean-cite-normalize': 'cite_normalize'
}


intg_settings = [
    intg_stvea_settings
]

intg_settings_keys = {
    'intg-STvEA': intg_stvea_settings_keys,
}
