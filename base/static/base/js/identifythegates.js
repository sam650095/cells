document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_results = await fetchAPI(
    "/api/preloadidentifythegates",
    0,
    csrftoken
  );
  console.log(preload_results);
  showdropdown(preload_results.adata_list);
});
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
  toggleLoading(true);
  const identifythegates_result = await fetchAPI(
    "/api/identifythegates",
    0,
    csrftoken
  );
  console.log(identifythegates_result);
  createTable(identifythegates_result.gate_df);
  toggleLoading(false);

  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("nextbtn").classList.remove("hidden");
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

  data.forEach((row, index) => {
    const tr = document.createElement("tr");
    tr.className = "bg-white border-b dark:bg-gray-800 dark:border-gray-700";

    const tdIndex = document.createElement("td");
    tdIndex.textContent = index + 1;
    tdIndex.className = "px-6 py-4";
    tr.appendChild(tdIndex);

    Object.values(row).forEach((value) => {
      const td = document.createElement("td");
      td.textContent = value !== null ? value : "N/A";
      td.className = "px-6 py-4";
      tr.appendChild(td);
    });

    tbody.appendChild(tr);
  });
}
