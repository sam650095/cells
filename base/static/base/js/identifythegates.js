let origindata = [];
let editdata = [];
document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_results = await fetchAPI(
    "/api/preloadidentifythegates",
    0,
    csrftoken
  );
  console.log(preload_results);
  showdropdown(preload_results.data.adata_list);
  // grab step
  const grabstep_chosen_rslt = await grabsteps(
    `/getSteps/identifythegates/chosen/`
  );
  const grabstep_rslt = await grabsteps(`/getSteps/identifythegates/process/`);

  if (
    (grabstep_chosen_rslt.message != "notfound") &
    (grabstep_rslt.message != "notfound")
  ) {
    stepped = true;
    document.getElementById("watchonly").classList.remove("hidden");
    document.getElementById("processbutton").classList.remove("hidden");
    document.getElementById("btnchanges").classList.remove("hidden");
    document.getElementById("nextbtn").classList.remove("hidden");
    const datatype = JSON.parse(grabstep_rslt.input_values["gate_df"]);
    console.log(datatype);
    // datatype.forEach((item, index) => {
    //   console.log(`**${index}**: ${JSON.stringify(item).replace(/"/g, "'")}`);
    // });
    createTable(datatype);

    document.getElementById("selected").textContent =
      grabstep_chosen_rslt.input_values["chosen_adata"];
    banned_operations();
  }
});
// banned operation
function banned_operations() {
  const processBtn = document.getElementById("processbutton");
  processBtn.disabled = true;
  document.getElementById("phenotyping_data").disabled = true;
  document.getElementById("confirmbtn").disabled = true;
}
function showdropdown(adata_list) {
  const preloadul = document.getElementById("preloadul");

  preloadul.innerHTML = "";
  for (let i = 0; i < adata_list.length; i++) {
    const li = document.createElement("li");
    li.innerHTML = `
                  <a onclick="changeselect('${adata_list[i]}')"
                  class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2"
              >${adata_list[i]}</a>
          `;
    preloadul.appendChild(li);
  }
}
async function changeselect(newText) {
  document.getElementById("selected").textContent = newText;
  const csrftoken = getCookie("csrftoken");
  const form = document.getElementById("identifythegatesform");
  const formData = new FormData(form);
  formData.append(
    "chosen_adata",
    document.getElementById("selected").textContent
  );
  const choseadata = await fetchAPI("/api/choseadata", formData, csrftoken);
  document.getElementById("chosen_results").textContent = choseadata.result;
  document.getElementById("processbutton").classList.remove("hidden");
}
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "processbutton");
  const identifythegates_result = await fetchAPI(
    "/api/identifythegates",
    0,
    csrftoken
  );
  createTable(identifythegates_result.data.gate_df);
  toggleLoading(false, "processbutton");

  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("btnchanges").classList.remove("hidden");
}
function createTable(data) {
  const table = document.getElementById("nametable");
  const thead = table.querySelector("thead tr");
  const tbody = document.getElementById("nametbody");

  thead.innerHTML = '<th scope="col" class="px-6 py-3">#</th>';
  tbody.innerHTML = "";

  Object.keys(data[0]).forEach((key) => {
    const th = document.createElement("th");
    th.textContent = key;
    th.scope = "col";
    th.className = "px-6 py-3";
    thead.appendChild(th);
  });
  origindata = JSON.parse(JSON.stringify(data));
  editdata = JSON.parse(JSON.stringify(data));
  data.forEach((row, index) => {
    const tr = document.createElement("tr");
    tr.className = "bg-white border-b hover:bg-gray-50 cursor-pointer";

    const tdIndex = document.createElement("td");
    tdIndex.textContent = index + 1;
    tdIndex.className = "px-6 py-4";
    tr.appendChild(tdIndex);

    Object.entries(row).forEach(([key, value]) => {
      const td = document.createElement("td");
      td.textContent = value !== null ? value : "N/A";
      td.className = "w-1/3 px-6 py-4";
      if (key != "Marker") td.onclick = () => handleClick(td, index, key);
      tr.appendChild(td);
    });

    tbody.appendChild(tr);
  });
}
function handleClick(td, rowIndex, columnName) {
  console.log("handleclick");
  if (td.querySelector("input")) {
    return;
  }
  const div = document.createElement("div");
  div.className = "flex justify-between";
  const currentValue = td.textContent;
  td.textContent = "";
  const input = document.createElement("input");

  input.type = "number";
  input.value = currentValue === "N/A" ? "" : currentValue;
  input.className =
    "disabled:bg-gray-300 disabled:cursor-not-allowed my-2 w-3/4 text-sm text-gray-900 border border-gray-300 rounded-lg bg-gray-50 focus:ring-sky-900 focus:border-sky-900";

  div.appendChild(input);
  const btn = document.createElement("button");
  btn.textContent = "V";

  btn.className =
    "disabled:bg-sky-700/50 disabled:border-sky-700/50 disabled:cursor-not-allowed disabled:text-white text-cyan-700 border border-sky-700 rounded-md hover:bg-sky-700 hover:text-white py-2 px-4 text-white text-sm font-semibold my-3";
  btn.onclick = function (event) {
    event.stopPropagation();
    saveEdit(rowIndex, columnName, input.value);
    td.innerHTML = "";
    td.textContent = input.value === "" ? "N/A" : input.value;
  };
  div.appendChild(btn);
  td.appendChild(div);
  input.focus();
  if (stepped) {
    input.disabled = true;
  }
}
function saveEdit(rowIndex, columnName, newValue) {
  editdata[rowIndex][columnName] = newValue;
  console.log("Updated editdata:", editdata);
}
async function confirmbtn() {
  const csrftoken = getCookie("csrftoken");
  const formData = new FormData();
  formData.append("editdata", JSON.stringify(editdata));
  const addvalueresult = await fetchAPI("/api/addvalue", formData, csrftoken);
  console.log(addvalueresult);
  document.getElementById("nextbtn").classList.remove("hidden");
}
