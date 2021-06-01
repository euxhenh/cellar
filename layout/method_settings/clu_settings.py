import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html


clu_leiden_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Leiden Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Resolution", html_for="slider"),
                        dcc.Slider(
                            id="clu-Leiden-resolution",
                            min=0.01, max=10, step=0.01,
                            value=1,
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("leidenalg.readthedocs.io/en/stable/reference.html#reference",
                               href="https://leidenalg.readthedocs.io/en/stable/reference.html#reference",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="clu-Leiden-settings",
    target="clu-Leiden-btn",
    trigger="click"
)

clu_leiden_settings_keys = {
    'clu-Leiden-resolution': 'resolution_parameter'
}


clu_settings = [clu_leiden_settings]
clu_settings_keys = {
    'clu-Leiden': clu_leiden_settings_keys,
}
